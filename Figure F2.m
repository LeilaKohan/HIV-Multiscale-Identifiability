clear all; close all; clc;


year_of_infection = [0,1,2,3,4,5,9,15,16,18,19,25,26,27,28,29,30,31,32,33,34,35,36,37,38]';

Estimated_HIV_incident_numbers = [20.77,65.39,65.13,131.07,131.37,85.37,85.15,49.53,49.89,58.6,59.19,56.40,56.49,53.70,42.21,41.10,40.28,38.89,38.07,38.66,38.95,38.68,37.85,37.30,36.19]';

year_of_infection_post = [39,40,41]'; Estimated_HIV_incident_numbers_post = [34.2,32.7,31.8]';

years_diagnose = (27:38)';

HIV_diagnose_cases = [46703,44324,42775,40945,40183,38931,39653,39543,39163,38048,37132,36350]'/1000;

years_diagnose_post = (39:43)';

HIV_diagnose_cases_post = [30.403,35.763,37.660,38.793,34.788]';

% Fitted parameters we found  
Lambda = 8.38663746616526;
beta = 0.306293687587148;
gamma = 0.120967014736454;
rho0 = 0.200000005673198;  
beta_d0 = 0.00119503317441228; % fitted tilde beta_d up to 2019

params=[Lambda,beta,gamma,rho0,beta_d0];

base_year=1981; 
dt=0.01; T=49; 
t1=2019;
t2=2030;  % simulate to 2030
goalI25=9.3;
goalI30=3.0; 
goalD25=9.588; 
goalD30=3.0;   % (in thousands)

% Plot color (consistent blue)
cBlue=[0 0.4470 0.7410];

% Rho-only optimization
obj_rho=@(rho_post) ...
    (get_incidence_at_year(2025,params,dt,T,base_year,rho_post,t1,t2)-goalI25)^2 + ...
    (get_incidence_at_year(2030,params,dt,T,base_year,rho_post,t1,t2)-goalI30)^2 + ...
    50*(get_diagnosis_at_year(2025,params,dt,T,base_year,rho_post,t1,t2)-goalD25)^2 + ...
    50*(get_diagnosis_at_year(2030,params,dt,T,base_year,rho_post,t1,t2)-goalD30)^2;

rho_post0=rho0+0.05; lb=rho0+0.01; ub=2.0;

opts=optimoptions('fmincon','Display','off','MaxFunctionEvaluations',3e4);

[rho_post_opt,~,~]=fmincon(obj_rho,rho_post0,[],[],[],[],lb,ub,[],opts);

[yrs,Inc,Diag]=project_full(params,dt,T,base_year,rho_post_opt,t1,t2);

i25=find(abs(yrs-2025)==min(abs(yrs-2025)),1); i30=find(abs(yrs-2030)==min(abs(yrs-2030)),1);

hit_rho = all([Inc(i25)<=goalI25, Inc(i30)<=goalI30, Diag(i25)<=goalD25, Diag(i30)<=goalD30]);

% joint optimization (rho + beta_d)
beta_d_post_opt = beta_d0;
if ~hit_rho
    obj_joint=@(x) ...
      (get_inc_beta(2025,params,dt,T,base_year,beta_d0,x(2),t1,t2,x(1))-goalI25).^2 + ...
      (get_inc_beta(2030,params,dt,T,base_year,beta_d0,x(2),t1,t2,x(1))-goalI30).^2 + ...
      50*(get_diag_beta(2025,params,dt,T,base_year,beta_d0,x(2),t1,t2,x(1))-goalD25).^2 + ...
      50*(get_diag_beta(2030,params,dt,T,base_year,beta_d0,x(2),t1,t2,x(1))-goalD30).^2;

    x0=[rho0+0.05; 0.7*beta_d0]; lb2=[rho0+0.01; 0]; ub2=[2.0; beta_d0];

    [x_opt,~,~]=fmincon(obj_joint,x0,[],[],[],[],lb2,ub2,[],opts);

    rho_post_opt=x_opt(1); beta_d_post_opt=x_opt(2);

    [yrs,Inc,Diag]=project_full_beta(params,dt,T,base_year,beta_d0,beta_d_post_opt,t1,t2,rho_post_opt);

    i25=find(abs(yrs-2025)==min(abs(yrs-2025)),1); i30=find(abs(yrs-2030)==min(abs(yrs-2030)),1);

end


% Incidence figure 
figure('Position',[100,100,1000,600]); hold on;

% 2019 change-point (hidden from legend)
xline(2019,'--k','LineWidth',1.2,'HandleVisibility','off');

% Model (blue), CDC <=2019 (purple), CDC 2020+ (red)
p1 = plot(yrs, Inc, '-', 'Color', cBlue, 'LineWidth', 2);
p2 = scatter(year_of_infection + base_year, Estimated_HIV_incident_numbers, ...
             80, [0.8 0 0.8], 'filled');
p3 = scatter(year_of_infection_post + base_year, Estimated_HIV_incident_numbers_post, ...
             80, [0.85 0 0], 'filled');

% Goal markers (not in legend)
scatter(2025, goalI25, 100, [0 0.6 0], 'v', 'filled', 'HandleVisibility','off');
scatter(2030, goalI30, 100, [0 0 0.8], 'v', 'filled', 'HandleVisibility','off');

% Legend with handles
legend([p1 p2 p3], ...
   {'Model predicted HIV incidences', ...
    'CDC reported incidences (\leq2019)', ...
    'CDC reported incidences (2020+)'}, ...
   'Location','northwest','AutoUpdate','off');

% Axes/labels/title
xlabel('Years of HIV Epidemic in USA','FontWeight','bold','FontSize',12);
ylabel('HIV Cases Per Year (in 1000s)','FontWeight','bold','FontSize',12);
title('HIV Incidences Projected to 2030');
set(gca,'FontSize',12,'LineWidth',1.5);
box on; grid on; xlim([1981,2032]);

% EHE text label
yl = ylim;
text(2013+0.2, yl(2)*0.75, 'Ending the HIV Epidemic Initiated', ...
     'FontSize',14,'Color',[0 0.5 0]);

% Goal annotations
arrow_offset = 40;
ax = gca;
ax_pos = get(ax, 'Position');
xlims = xlim(ax);
ylims = ylim(ax);
x2fig = @(x) ax_pos(1) + (x - xlims(1)) / diff(xlims) * ax_pos(3);
y2fig = @(y) ax_pos(2) + (y - ylims(1)) / diff(ylims) * ax_pos(4);

annotation('textarrow', ...
    [x2fig(2025 - 2.5), x2fig(2025)], ...
    [y2fig(goalI25 + arrow_offset), y2fig(goalI25)], ...
    'String','2025 Year Target: 9.3k', ...
    'FontSize',12,'Color',[0 0.5 0]);

annotation('textarrow', ...
    [x2fig(2030 - 3), x2fig(2030)], ...
    [y2fig(goalI30 + arrow_offset), y2fig(goalI30)], ...
    'String','2030 Year Target: 3k', ...
    'FontSize',12,'Color',[0 0 0.6]);



% Diagnoses figure 
figure('Position',[100,100,1000,600]); hold on;

% 2019 change-point (hidden from legend)
xline(2019,'--k','LineWidth',1.2,'HandleVisibility','off');

% Model (blue), CDC <=2019 (purple), CDC 2020+ (red)
p1 = plot(yrs, Diag, '-', 'Color', cBlue, 'LineWidth', 2);
p2 = scatter(years_diagnose + base_year, HIV_diagnose_cases, ...
             80, [0.8 0 0.8], 'filled');
p3 = scatter(years_diagnose_post + base_year, HIV_diagnose_cases_post, ...
             80, [0.85 0 0], 'filled');

% Goal markers (not in legend)
scatter(2025, goalD25, 100, [0 0.6 0], 'v', 'filled', 'HandleVisibility','off');
scatter(2030, goalD30, 100, [0 0 0.8], 'v', 'filled', 'HandleVisibility','off');

% Legend with handles
legend([p1 p2 p3], ...
   {'Model predicted HIV diagnoses', ...
    'CDC reported diagnoses (≤2019)', ...
    'CDC reported diagnoses (2020+)'}, ...
   'Location','northwest','AutoUpdate','off');

% Axes/labels/title
xlabel('Years of HIV Epidemic in USA','FontWeight','bold','FontSize',12);
ylabel('HIV Diagnoses Per Year (in 1000s)','FontWeight','bold','FontSize',12);
title('HIV Diagnoses Projected to 2030');
set(gca,'FontSize',12,'LineWidth',1.5);
box on; grid on; xlim([1981,2032]);

yl = ylim;
text(2013 + 0.2, yl(2) * 0.75, ...
'Ending the HIV Epidemic Initiated', ...
'FontSize', 14, 'Color', [0 0.5 0]);
% Add triangles
scatter(2025, goalD25, 100, [0 0.6 0], 'v', 'filled');
scatter(2030, goalD30, 100, [0 0 0.8], 'v', 'filled');
% arrow annotations
ax = gca;
ax_pos = get(ax, 'Position');
xlim_ax = xlim(ax); ylim_ax = ylim(ax);
x2fig = @(x) ax_pos(1) + (x - xlim_ax(1)) / diff(xlim_ax) * ax_pos(3);
y2fig = @(y) ax_pos(2) + (y - ylim_ax(1)) / diff(ylim_ax) * ax_pos(4);
arrow_offset = 40;
% 2025 Arrow
annotation('textarrow', ...
[x2fig(2025 - 2.5), x2fig(2025)], ...
[y2fig(goalD25 + arrow_offset), y2fig(goalD25)], ...
'String', '2025 Year Target: 9.6k', ...
'FontSize', 12, 'Color', [0 0.5 0]);
% 2030 Arrow
annotation('textarrow', ...
[x2fig(2030 - 3), x2fig(2030)], ...
[y2fig(goalD30 + arrow_offset), y2fig(goalD30)], ...
'String', '2030 Year Target: 3k', ...
'FontSize', 12, 'Color', [0 0 0.6]);
hold off;

%  One graph: rho(t) and \tilde{\beta}_d(t) 
tt = 1981:2030;

rho_vals = arrayfun(@(t) piecewise_linear_rho(t, rho0, rho_post_opt, t1, t2), tt);

% If beta_d_post_opt exists use it; otherwise show a flat line at beta_d0
if exist('beta_d_post_opt','var')
    beta_d_vals = arrayfun(@(t) piecewise_linear_rho(t, beta_d0, beta_d_post_opt, t1, t2), tt);
else
    beta_d_vals = beta_d0 * ones(size(tt));
end



% rho(t) and \tilde{\beta}_d(t)
tt = 1981:2030;
cBlue = [0 0.4470 0.7410];

% series
rho_vals   = arrayfun(@(t) piecewise_linear_rho(t, rho0, rho_post_opt, t1, t2), tt);
if exist('beta_d_post_opt','var')
    beta_d_vals = arrayfun(@(t) piecewise_linear_rho(t, beta_d0, beta_d_post_opt, t1, t2), tt);
else
    beta_d_vals = beta_d0*ones(size(tt));
end

% helpful points (2019/2030)
rho_2019 = piecewise_linear_rho(2019, rho0, rho_post_opt, t1, t2);
rho_2030 = piecewise_linear_rho(2030, rho0, rho_post_opt, t1, t2);
bd_2019  = piecewise_linear_rho(2019, beta_d0, beta_d_vals(end), t1, t2);
bd_2030  = piecewise_linear_rho(2030, beta_d0, beta_d_vals(end), t1, t2);

figure('Position',[100,100,1200,480]); hold on;

% ---- left axis: rho(t)
yyaxis left
h1 = plot(tt, rho_vals, '-', 'LineWidth', 2.4, 'Color', cBlue);
padL = max(1e-3, 0.05*range(rho_vals));
ylim([min(rho_vals)-padL, max(rho_vals)+padL]);
ylabel('$\rho(t)\ \mathrm{(yr^{-1})}$', 'Interpreter','latex','FontSize',16);
set(gca,'YColor','k')

% markers (rho)
plot([2019 2030],[rho_2019 rho_2030], 'o', 'MarkerSize',6, ...
     'MarkerFaceColor', cBlue, 'MarkerEdgeColor', cBlue, 'LineStyle','none');

% ---- right axis: \tilde{\beta}_d(t) in x10^{-3}
yyaxis right
h2 = plot(tt, 1e3*beta_d_vals, '--', 'LineWidth', 2.4, 'Color', cBlue);
padR = max(1e-3, 0.08*range(1e3*beta_d_vals));
ylim([max(0,min(1e3*beta_d_vals)-padR), max(1e3*beta_d_vals)+padR]);
ylabel('$\tilde{\beta}_d(t)\ (\times 10^{-3})$', 'Interpreter','latex','FontSize',16);
set(gca,'YColor','k')



% markers (beta_d)
plot([2019 2030], 1e3*[bd_2019 bd_2030], '^', 'MarkerSize',6, ...
     'MarkerFaceColor', cBlue, 'MarkerEdgeColor', cBlue, 'LineStyle','none');

% common cosmetics
xlim([1981 2032])
grid on; box on; set(gca,'LineWidth',1.2,'FontSize',12,'TickDir','out')
xlabel('Year', 'FontWeight','bold')
title({'Diagnosis rate and diagnosed-transmission'}, 'FontWeight','bold')

% vertical line + shaded transition window
ylL = ylim; yyaxis left; ylL = ylim;      % lock left limits
xline(2019,'--','Color',[0.25 0.25 0.25],'LineWidth',1.2,'HandleVisibility','off')
% subtle shading for 2019–2030
yl = ylim;
patch([2019 2030 2030 2019],[yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
      'FaceAlpha',0.07,'EdgeColor','none','HandleVisibility','off');

% legend (LaTeX)
yyaxis left
lgd = legend([h1 h2], {'$\rho(t)$ (left axis)', '$\tilde{\beta}_d(t)$ (right axis)'}, ...
             'Interpreter','latex','Location','northwest','FontSize', 16);
set(lgd,'Box','off')

% optional: save hi-res
% exportgraphics(gcf,'rho_beta_single.png','Resolution',300);



% R0(t) in same blue
[R0_years,R0_series,~,~,~]=compute_R0_series(params,dt,T,base_year,beta_d0,beta_d_post_opt,t1,t2);

figure('Position',[100,100,950,420]);
hold on;
plot(R0_years(:),R0_series(:),'-','Color',cBlue,'LineWidth',2.2);

xline(2019,'--','Color',[0 0 0],'LineWidth',1.2); yline(1,'--','Color',[0 0.6 0],'LineWidth',1.2);
xlabel('Year','FontWeight','bold'); 
ylabel('Basic reproduction number R_0(t)','FontSize',14);
title({'Refitting $\tilde{\beta}_d$ and $\rho$ after 2019'}, ...
      'Interpreter','latex','FontSize',18);
grid on;
box on;

% Helpers 
function inc=get_incidence_at_year(y,params,dt,T,base,rho_post,t1,t2)
    [yrs,inc,~]=project_full(params,dt,T,base,rho_post,t1,t2); [~,i]=min(abs(yrs-y)); inc=inc(i);
end

function diag=get_diagnosis_at_year(y,params,dt,T,base,rho_post,t1,t2)
    [yrs,~,diag]=project_full(params,dt,T,base,rho_post,t1,t2); [~,i]=min(abs(yrs-y)); diag=diag(i);
end

function inc=get_inc_beta(y,params,dt,T,base,beta_d0,beta_d_post,t1,t2,rho_post)
    [yrs,inc,~]=project_full_beta(params,dt,T,base,beta_d0,beta_d_post,t1,t2,rho_post); [~,i]=min(abs(yrs-y)); inc=inc(i);
end

function diag=get_diag_beta(y,params,dt,T,base,beta_d0,beta_d_post,t1,t2,rho_post)
    [yrs,~,diag]=project_full_beta(params,dt,T,base,beta_d0,beta_d_post,t1,t2,rho_post); [~,i]=min(abs(yrs-y)); diag=diag(i);
end

function [years,HIV_Inc,HIV_Diag]=project_full(params,dt,T_end,base_year,rho_post,t_start,t_end)

    Lambda=params(1); 
    beta=params(2); 
    gamma=params(3); rho0=params(4); 
    beta_d0=params(5);
    d_tau=0.52;
    kappa=52; 
    Treat_Age=20*52;
    tforward=0:dt:T_end; 
    treat_time=0:d_tau:Treat_Age;
    fp=[33.84,0.0087,0.00036,0.085,0.45,0.96];

    dw=@(t,y)[fp(1)-fp(2)*y(1)*y(3)-fp(3)*y(1); fp(2)*y(1)*y(3)-fp(4)*y(2); fp(5)*y(2)-fp(6)*y(3)];
    [~,ys]=ode23s(dw,treat_time,[167;85;69643]); Vvec=ys(:,3);

    N=length(tforward); M=length(treat_time);

    S=zeros(N,1); 
    I=zeros(N,1); 
    A=zeros(N,1); 
    D=zeros(N,1); 
    d_old=zeros(M,1);
    d_new=zeros(M,1);
    S(1)=820; 
    I(1)=42; 
    A(1)=0.337;
    D(1)=0;
    Bd=zeros(N,1); 
    Nt=zeros(N,1); 
    Nt(1)=S(1)+I(1)+D(1);

    for k=1:N-1
        tcurr=tforward(k)+base_year;

        rhot=piecewise_linear_rho(tcurr,rho0,rho_post,t_start,t_end);

        beta_d_vec=beta_d0*Vvec;

        Nt(k)=S(k)+I(k)+D(k); 

        if Nt(k)==0, Nt(k)=1;
        end

        Bd(k)=d_tau*trapz(beta_d_vec.*d_old)/Nt(k);

        S(k+1)=(S(k)+dt*Lambda)/(1+dt*(beta*I(k)/Nt(k)+0.02));

        I(k+1)=(I(k)+dt*(beta*S(k+1)*I(k)/Nt(k)+S(k+1)*Bd(k)))/(1+(0.02+0.002+rhot)*dt);

        d_new(1)=(1/kappa)*rhot*I(k+1); d_new(2:end)=d_old(1:end-1)./(1+dt*0.02);

        A(k+1)=(A(k)+dt*(gamma*I(k+1)))/(1+dt*(0.02+0.002));

        D(k+1)=trapz(d_new)*d_tau; d_old=d_new;
    end

    HIV_Inc=beta*S(1:end-1).*I(1:end-1)./Nt(1:end-1)+S(1:end-1).*Bd(1:end-1);

    years=tforward(1:end-1)+base_year;

    rho_vec=arrayfun(@(tt) piecewise_linear_rho(tt,rho0,rho_post,t_start,t_end), years);

    HIV_Diag=rho_vec(:).*I(1:end-1);
end

function [years,HIV_Inc,HIV_Diag]=project_full_beta(params,dt,T_end,base_year,beta_d0,beta_d_post,t_start,t_end,rho_post)

    Lambda=params(1);
    beta=params(2); 
    gamma=params(3); 
    rho0=params(4);

    d_tau=0.52; 
    kappa=52; 
    Treat_Age=20*52;
    tforward=0:dt:T_end;
    treat_time=0:d_tau:Treat_Age;
   
    fp=[33.84,0.0087,0.00036,0.085,0.45,0.96];

    dw=@(t,y)[fp(1)-fp(2)*y(1)*y(3)-fp(3)*y(1); fp(2)*y(1)*y(3)-fp(4)*y(2); fp(5)*y(2)-fp(6)*y(3)];

    [~,ys]=ode23s(dw,treat_time,[167;85;69643]); Vvec=ys(:,3);

    N=length(tforward); M=length(treat_time);

    S=zeros(N,1); 
    I=zeros(N,1); 
    A=zeros(N,1);
    D=zeros(N,1); 
    d_old=zeros(M,1); 
    d_new=zeros(M,1);

    S(1)=820;
    I(1)=42; 
    A(1)=0.337;
    D(1)=0; 

    Bd=zeros(N,1);
    Nt=zeros(N,1); 
    Nt(1)=S(1)+I(1)+D(1);

    beta_d_fun=@(t) piecewise_linear_rho(t,beta_d0,beta_d_post,t_start,t_end);

    for k=1:N-1
        tcurr=tforward(k)+base_year;
        rhot=piecewise_linear_rho(tcurr,rho0,rho_post,t_start,t_end);
        beta_d_tilde=beta_d_fun(tcurr); beta_d_vec=beta_d_tilde*Vvec;
        Nt(k)=S(k)+I(k)+D(k);

        if Nt(k)==0, Nt(k)=1;
        end
        Bd(k)=d_tau*trapz(beta_d_vec.*d_old)/Nt(k);
        S(k+1)=(S(k)+dt*Lambda)/(1+dt*(beta*I(k)/Nt(k)+0.02));
        I(k+1)=(I(k)+dt*(beta*S(k+1)*I(k)/Nt(k)+S(k+1)*Bd(k)))/(1+(0.02+0.002+rhot)*dt);
        d_new(1)=(1/kappa)*rhot*I(k+1); d_new(2:end)=d_old(1:end-1)./(1+dt*0.02);
        A(k+1)=(A(k)+dt*(gamma*I(k+1)))/(1+dt*(0.02+0.002));
        D(k+1)=trapz(d_new)*d_tau; d_old=d_new;
    end

    HIV_Inc=beta*S(1:end-1).*I(1:end-1)./Nt(1:end-1)+S(1:end-1).*Bd(1:end-1);
    years=tforward(1:end-1)+base_year;
    rho_vec=arrayfun(@(tt) piecewise_linear_rho(tt,rho0,rho_post,t_start,t_end), years);
    HIV_Diag=rho_vec(:).*I(1:end-1);
end

function val=piecewise_linear_rho(t,rho0,rho_post,t1,t2)
    if t<=t1, val=rho0;
    elseif t<=t2, val=rho0+(rho_post-rho0)*(t-t1)/(t2-t1);
    else, val=rho_post;
    end
end

function [years,R0_series,R0_base,R0_coef,R0_pre2019]=compute_R0_series(params,dt,T_end,base_year,beta_d0,beta_d_post,t_start,t_end)

    beta=params(2); gamma=params(3); rho=params(4); mu=0.02; kappa=52; d_tau=0.52; Treat_Age=20*52;

    t_tau=(0:d_tau:Treat_Age)'; fp=[33.84,0.0087,0.00036,0.085,0.45,0.96];

    init=[167;85;69643];

    dydt=@(t,y)[fp(1)-fp(2)*y(1)*y(3)-fp(3)*y(1); fp(2)*y(1)*y(3)-fp(4)*y(2); fp(5)*y(2)-fp(6)*y(3)];

    [~,ys]=ode23s(dydt,t_tau,init); Vtau=ys(:,3);

    C=trapz(t_tau, Vtau .* exp(-(mu/kappa).*t_tau));

    denom=(mu+gamma+rho); 
    R0_base=beta/denom;
    R0_coef=(rho/(kappa*denom))*C; 
    R0_pre2019=R0_base+R0_coef*beta_d0;
    tforward=(0:dt:T_end)';
    years=tforward+base_year; 
    R0_series=nan(size(years));

    for i=1:numel(years)
        t=years(i);
        beta_d_t = piecewise_linear_rho(t,beta_d0,beta_d_post,t_start,t_end);
        R0_series(i)=R0_base + R0_coef*beta_d_t;
    end
end
