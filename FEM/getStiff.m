function [K, d_K] = getStiff(domain)
%Compute global stiffness matrix K

%row and column indices
neq = size(domain.nodes, 1) - sum(domain.nodes);
r = zeros(3*neq - 2, 1);
c = zeros(3*neq - 2, 1);
Kv = zeros(3*neq - 2, 1);
d_K = zeros(neq, neq, domain.N_el);
for i = 1:(domain.N_el + 1)
    
   if ~domain.nodes(i)
       
       if ~domain.nodes(1)
           
          if i == 1
             r(1:2) = i;
             c(1:2) = [1; 2];
             Kv(1:2) = (1/domain.l)*[domain.conductivity(1); -domain.conductivity(1)];
             d_K(1, 1:2, 1) = (1/domain.l)*[1, -1];
             
          elseif i == domain.N_el + 1 || (i == domain.N_el && domain.nodes(end))
              %last natural node
              r((3*neq - 3):(3*neq - 2)) = neq;
              c((3*neq - 3):(3*neq - 2)) = [neq - 1; neq];
              if i == domain.N_el
                  %last node is essential
                  Kv((3*neq - 3):(3*neq - 2)) = (1/domain.l)*[-domain.conductivity(i - 1);...
                      domain.conductivity(i - 1) + domain.conductivity(i)];
                  d_K(neq, (neq - 1):neq, i - 1) = (1/domain.l)*[-1, 1];
                  d_K(neq, neq, i) = (1/domain.l);
              else
                  %last node is natural
                  Kv((3*neq - 3):(3*neq - 2)) = (1/domain.l)*[-domain.conductivity(i - 1);...
                      domain.conductivity(i - 1)];
                  d_K(neq, (neq - 1):neq, i - 1) = (1/domain.l)*[-1, 1];
              end
          else
              r(3*(i - 1):(3*i - 1)) = i;
              c(3*(i - 1):(3*i - 1)) = [i - 1; i; i + 1];
              Kv(3*(i - 1):(3*i - 1)) = (1/domain.l)*[-domain.conductivity(i - 1);...
                  domain.conductivity(i - 1) + domain.conductivity(i); - domain.conductivity(i)];
              d_K(i, (i - 1):i, i - 1) = (1/domain.l)*[-1, 1];
              d_K(i, i:(i + 1), i) = (1/domain.l)*[1, -1];
          end  
          
       else
           %first node is essential
           if i == 1
               
               continue     %skip first node
               
           elseif i == 2
             r(1:2) = i - 1;
             c(1:2) = [1; 2];
             Kv(1:2) = (1/domain.l)*[domain.conductivity(i - 1) + domain.conductivity(i); -domain.conductivity(i)];
             d_K(i - 1, 1, i - 1) = (1/domain.l);
             d_K(i - 1, 1:2, i) = (1/domain.l)*[1, -1];
             
          elseif i == domain.N_el + 1 || (i == domain.N_el && domain.nodes(end))
              %last natural node
              r((3*neq - 3):(3*neq - 2)) = neq;
              c((3*neq - 3):(3*neq - 2)) = [neq - 1; neq];
              if i == domain.N_el
                  %last node is essential
                  Kv((3*neq - 3):(3*neq - 2)) = (1/domain.l)*[-domain.conductivity(i - 1);...
                      domain.conductivity(i - 1) + domain.conductivity(i)];
                  d_K(neq, (neq - 1):neq, i - 1) = (1/domain.l)*[-1, 1];
                  d_K(neq, neq, i) = (1/domain.l);
              else
                  %last node is natural
                  Kv((3*neq - 3):(3*neq - 2)) = (1/domain.l)*[-domain.conductivity(i - 1);...
                      domain.conductivity(i - 1)];
                  d_K(neq, (neq - 1):neq, i - 1) = (1/domain.l)*[-1, 1];
              end
          else
              r(3*(i - 2):(3*(i - 1) - 1)) = i - 1;
              c(3*(i - 2):(3*(i - 1) - 1)) = [i - 2; i - 1; i];
              Kv(3*(i - 2):(3*(i - 1) - 1)) = (1/domain.l)*[-domain.conductivity(i - 1);...
                  domain.conductivity(i - 1) + domain.conductivity(i); - domain.conductivity(i)];
              d_K(i - 1, (i - 2):(i - 1), i - 1) = (1/domain.l)*[-1, 1];
              d_K(i - 1, (i - 1):i, i) = (1/domain.l)*[1, -1];
          end 
           
       end
       
   end
    
end
K = sparse(r, c, Kv);
    
end

