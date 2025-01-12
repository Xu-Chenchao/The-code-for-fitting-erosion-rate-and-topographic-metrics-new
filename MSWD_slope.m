%% calculate MSWD ��ɢ��Լ���ĺ������ߺ����ݵ㱾��������
function MSWD = MSWD(E_pred,topo_metrics_pred,E,E_error,topo_metrics,topo_metrics_error,data5)

% E_pred and topo_metrics_pred stand for the variables x and y
% E and topo_metrics stand for the observed point. E_error and topo_metrics_error
% are the error of the E and topo_metrics
% data5 is all data (including all bins)

n = length(E); % �в����
d = 2; % ���ɶ�����Ϊ2
r2w = []; %���ڴ洢��Ȩ�صĲв�
% emax = mean(data5(:,1)) + 10*std(data5(:,1));%�ù۲���ʴ���ʵ�ƽ��ֵ����10���ı�׼����ΪSc���ϵ�erosion rate��Ԥ��ֵ
% emax = (mean(data5(:,1)) + 4*std(data5(:,1)))/1000;%�ù۲���ʴ���ʵ�ƽ��ֵ����10���ı�׼����ΪSc���ϵ�erosion rate��Ԥ��ֵ m/yr
emax = (mean(E) + 4*std(E))/1000;%�ù۲���ʴ���ʵ�ƽ��ֵ����10���ı�׼����ΪSc���ϵ�erosion rate��Ԥ��ֵ m/yr �÷���bin�е���ʴ����ƽ��ֵ
for k = 1:length(E)
    % �������������ʵ������֮��Ĳв�
    if topo_metrics(k) < topo_metrics_pred(end-40)
        residuals = interp1(topo_metrics_pred,E_pred,topo_metrics(k)) - E(k);  %����MSWD�Ĳв�  erosion rate ��MSWD
%         interp_func = scatteredInterpolant(topo_metrics_pred,E_pred');
%         E_value = interp_func(topo_metrics(k));
%         residuals = E_value - E(k);
    elseif topo_metrics(k) > topo_metrics_pred(end-40)
%     TF = isnan(residuals);
%     if TF == 1
%         residuals = emax - E(k);
          residuals = 0;
%     end
    end
   % residuals = interp1(E_pred,topo_metrics_pred,E(k)) - topo_metrics(k);  %����MSWD�Ĳв�  erosion rate ��MSWD
%     slo = topo_metrics_error(k)/E_error(k);
%     x = E_pred;
%     y = E - slo*(x-topo_metrics);  %���۲���ֱ�ߵı��ʽ��x��Χ���ĺ�E_predһ��
%     a = [];
%     for i = 1:length(E)
%         res = (topo_metrics_pred(i) - y(i))^2;
%         a(i) = res;
%     end
%     idx = find(min(a));
%     % ��������
%     x_intersect = E_pred(idx); % ��ȡ����� x ����
%     y_intersect = topo_metrics_pred(idx); % ��ȡ����� y ����
%     % ���������ڽ��㴦������б�ʣ������㴦��б�ʣ�
%     if idx == 1
%        slope_at_intersection = (topo_metrics_pred(idx+1) - topo_metrics_pred(idx)) / (E_pred(idx+1) - E_pred(idx));
%     elseif idx == length(x_fit)
%        slope_at_intersection = (topo_metrics_pred(idx) - topo_metrics_pred(idx-1)) / (E_pred(idx) - E_pred(idx-1));
%     else
%        slope_at_intersection = (topo_metrics_pred(idx+1) - topo_metrics_pred(idx-1)) / (E_pred(idx+1) - E_pred(idx-1));
%     end
% 
%     % weights = ones(size(residuals)); % ����Ȩ�أ�����򵥵�����Ϊ1��
%     % weights = 1/((slope_at_intersection^2)*(E_error(k)^2)+(1^2)*(topo_metrics(k)^2));
%     weights = 1/((E_error(k)^2)+(slope_at_intersection^2)*(topo_metrics_error(k)^2)); %�ο�adams�Ĵ���д��
%     % ��Ȩ�صĲв�
%     if k == 1
%         weights = 1/(E_error(k)^2 + (E_pred(k+1)-E_pred(k))/(topo_metrics_pred(k+1)-topo_metrics_pred(k))*topo_metrics_error(k)^2);
%     elseif k == length(E)
%         weights = 1/(E_error(k)^2 + (E_pred(k)-E_pred(k-1))/(topo_metrics_pred(k)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
%     else
%         weights = 1/(E_error(k)^2 + (E_pred(k+1)-E_pred(k-1))/(topo_metrics_pred(k+1)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
%     end
% E unit m/Myr
    if k == 1
        weights = 1/((E_error(k)*1000000)^2 + (E_pred(k+1)-E_pred(k))*1000000/(topo_metrics_pred(k+1)-topo_metrics_pred(k))*topo_metrics_error(k)^2);
    elseif k == length(E)
        weights = 1/((E_error(k)*1000000)^2 + (E_pred(k)-E_pred(k-1))*1000000/(topo_metrics_pred(k)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
    else
        weights = 1/((E_error(k)*1000000)^2 + (E_pred(k+1)-E_pred(k-1))*1000000/(topo_metrics_pred(k+1)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
    end
    r2w(k) = ((residuals*10^6)^2)*weights;
end
% mean_squared_residuals = sum(r2w) / (n-d); % ����ƽ��ƽ���в� %�����NAN�����ֵ�ĺ�
mean_squared_residuals = sum(r2w) / (nnz(r2w)-d); % ����ƽ��ƽ���в� %�����NAN�����ֵ�ĺ�
% mean_squared_residuals = nansum(r2w) / (nnz(~isnan(r2w))-d); % ����ƽ��ƽ���в� %�����NAN�����ֵ�ĺ�
MSWD = mean_squared_residuals; % MSWD��Ϊƽ��ƽ���в�
end