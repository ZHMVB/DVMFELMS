function [f_output,f_middle_output] = f_Sz(f_X,f_middle,f_S1,f_S2,N1,N2)
%   非线性次级通路, f_output为输出，f_middle为中间变量
%   LNL模型
f_x1 = f_X(1:N1)' * f_S1;
f_x1 = 3.3*tanh(0.3*f_x1);
f_middle = [f_x1;f_middle(1:end-1)];
f_output = f_middle(1:N2)' * f_S2;
f_middle_output = f_middle;
end

