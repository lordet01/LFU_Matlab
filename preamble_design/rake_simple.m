close('all');

 A = sin(-pi:0.1:pi)';
 A = [A;zeros(100,1)];
 Td = 10;
 B = 0.5 * [zeros(Td,1); A(1:end-Td)];

 %Multi-path interference
 C = A + B;
 
 %Rake Receiver
 D = C + 0.5 * [C(Td:end);zeros(Td-1,1)];
 D = D ./ 2;
  
 plot(A); hold on; plot(B,'r'); plot(C,'k');  plot(D,'c'); hold off;
