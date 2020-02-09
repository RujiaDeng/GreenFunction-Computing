function rootBessel = GCL_rootBessel()    

% ���n�ױ��������������(0-9)
% nΪ�������������� 
% NΪҪ����������
% a��b��Ԫ��Ϊa-1�ױ����������ĵ�b�����
n = 9;
N = 50;
j = zeros(n+1, N);    % �����������ĸ�
incr = 4.0;
for v = 0 : n
   h = v + 1.9*v^(1/3)+1;
   if (v == 0)             % 0�ױ����������ĵ�һ�����
       j(v+1,1) = fzero(@(x)besselj(v,x),2);
   else                    % 1�׼����Ͻױ����������ĵ�һ�����
       j(v+1,1) = fzero(@(x)besselj(v,x),h);
   end
   for s = 2 : N           % �����������ĵ�2������������
       j(v+1,s) = fzero(@(x)besselj(v,x),j(v+1,s-1)+incr);
   end    
end

rootBessel = j;