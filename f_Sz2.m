function [f_output,f_middle_output] = f_Sz2(f_X,f_middle1,f_middle2,f_S1,f_S2,N1,N2)
%   ����������μ�ͨ·s_e, f_outputΪ�����f_middleΪ�м����
%   LNLģ��
f_x1 = f_X(1:N2)' * f_S2;
f_x1 = f_x1 * 3.3 * 0.3 * (1 - (tanh(0.3 * f_middle1(N2)))^2 );
f_middle2 = [f_x1;f_middle2(1:end-1)];
f_output = f_middle2(1:N1)' * f_S1;
f_middle_output = f_middle2;
end

