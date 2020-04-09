load Reluplex_P4.mat;
Reluplex_Res = Res;
Reluplex_VT = VT;
load P4_exact_star.mat;
exact_star_Res = Res;
exact_star_VT = VT; 
load P4_approx_star.mat;
approx_star_Res = Res; 
approx_star_VT = VT; 
load P4_approx_zono.mat;
approx_zono_Res = Res;
approx_zono_VT = VT;
load P4_abs_dom.mat;
abs_dom_Res = Res; 
abs_dom_VT = VT; 

exact_star_imp = zeros(5,9);
approx_star_imp = zeros(5,9);
approx_zono_imp = zeros(5,9);
abs_dom_imp = zeros(5,9);
for i=1:5
    for j=1:9
        exact_star_imp(i,j) = Reluplex_VT(i,j)/exact_star_VT(i,j);
        approx_star_imp(i,j) = Reluplex_VT(i,j)/approx_star_VT(i,j);
        approx_zono_imp(i,j) = Reluplex_VT(i,j)/approx_zono_VT(i,j);
        abs_dom_imp(i,j) = Reluplex_VT(i,j)/abs_dom_VT(i,j);
    end
end
exact_star_average_imp = sum(exact_star_imp, 'all')/45;
approx_star_average_imp = sum(approx_star_imp, 'all')/45;
approx_zono_average_imp = sum(approx_zono_imp, 'all')/45;
abs_dom_average_imp = sum(abs_dom_imp, 'all')/45;

Reluplex_safes = 0;
exact_star_safes = 0;
approx_star_safes = 0;
approx_zono_safes = 0;
abs_dom_safes = 0;

for i=1:5
    for j=1:9
        
        if strcmp(Reluplex_Res{i,j}, 'UNSAT')
            Reluplex_safes = Reluplex_safes + 1;
        end
        if strcmp(exact_star_Res{i,j}, 'UNSAT')
            exact_star_safes = exact_star_safes + 1;
        end
        if strcmp(approx_star_Res{i,j}, 'UNSAT')
            approx_star_safes = approx_star_safes + 1;
        end
        if strcmp(approx_zono_Res{i,j}, 'UNSAT')
            approx_zono_safes = approx_zono_safes + 1;
        end
        if strcmp(abs_dom_Res{i,j}, 'UNSAT')
            abs_dom_safes = abs_dom_safes + 1;
        end  
    end
end


% print to screen
fprintf('\n===============================================VERIFICATION RESULTS FOR PROPERTY P4====================================================\n');
fprintf('\n=======================================================================================================================================\n');
fprintf('\n%s%15s%25s%25s%25s%30s', 'Network ID', 'Reluplex','Exact-Star', 'Approx-star', 'Zonotope', 'Abstract-Domain');
fprintf('\n%18s%8s%12s%8s%8s%10s%8s%8s%10s%8s%8s%10s%8s%8s', 'Res', 'VT', 'Res', 'VT', 'Imp', 'Res', 'VT', 'Imp', 'Res', 'VT', 'Imp', 'Res', 'VT', 'Imp');
for i=1:5
    for j=1:9
        fprintf('\n%-14s%-8s%-7.1f%10s%7.1f%9.1f%10s%7.1f%10.1f%8s%9.1f%10.1f%7s%9.1f%10.1f',sprintf('N_%d%d', i, j), Reluplex_Res{i,j}, Reluplex_VT(i,j),exact_star_Res{i,j}, exact_star_VT(i,j),Reluplex_VT(i,j)/exact_star_VT(i,j), approx_star_Res{i,j}, approx_star_VT(i,j),Reluplex_VT(i,j)/approx_star_VT(i,j), approx_zono_Res{i,j}, approx_zono_VT(i,j),Reluplex_VT(i,j)/approx_zono_VT(i,j), abs_dom_Res{i,j}, abs_dom_VT(i,j),Reluplex_VT(i,j)/abs_dom_VT(i,j));
    end
end
% 
fprintf('\nTotal VT:%20.2f%20.2f%25.2f%26.2f%26.2f', sum(Reluplex_VT, 'all'), sum(exact_star_VT, 'all'), sum(approx_star_VT, 'all'), sum(approx_zono_VT, 'all'), sum(abs_dom_VT, 'all'));
fprintf('\nAverage Imp:%43.1f%28.1f%27.1f%26.1f', exact_star_average_imp, approx_star_average_imp, approx_zono_average_imp, abs_dom_average_imp);
fprintf('\nSafe Nets:%5d/%2d%18d/%2d%24d/%2d%22d/%2d%24d/%2d',Reluplex_safes, 45, exact_star_safes, 45, approx_star_safes, 45, approx_zono_safes, 45, abs_dom_safes, 45);



% generate txt file
fid = fopen('P4_tab.txt', 'wt');
fprintf(fid,'\n===============================================VERIFICATION RESULTS FOR PROPERTY P4====================================================\n');
fprintf(fid,'\n=======================================================================================================================================\n');
fprintf(fid,'\n%s%15s%25s%25s%25s%30s', 'Network ID', 'Reluplex','Exact-Star', 'Approx-star', 'Zonotope', 'Abstract-Domain');
fprintf(fid,'\n%18s%8s%12s%8s%8s%10s%8s%8s%10s%8s%8s%10s%8s%8s', 'Res', 'VT', 'Res', 'VT', 'Imp', 'Res', 'VT', 'Imp', 'Res', 'VT', 'Imp', 'Res', 'VT', 'Imp');

for i=1:5
    for j=1:9
        fprintf(fid,'\n%-14s%-8s%-7.1f%10s%7.1f%9.1f%10s%7.1f%10.1f%10s%7.1f%10.1f%10s%6.1f%10.1f',sprintf('N_%d%d', i, j), Reluplex_Res{i,j}, Reluplex_VT(i,j),exact_star_Res{i,j}, exact_star_VT(i,j),Reluplex_VT(i,j)/exact_star_VT(i,j), approx_star_Res{i,j}, approx_star_VT(i,j),Reluplex_VT(i,j)/approx_star_VT(i,j), approx_zono_Res{i,j}, approx_zono_VT(i,j),Reluplex_VT(i,j)/approx_zono_VT(i,j), abs_dom_Res{i,j}, abs_dom_VT(i,j),Reluplex_VT(i,j)/abs_dom_VT(i,j));
    end
end
% 
fprintf(fid,'\nTotal VT:%20.2f%20.2f%25.2f%26.2f%26.2f', sum(Reluplex_VT, 'all'), sum(exact_star_VT, 'all'), sum(approx_star_VT, 'all'), sum(approx_zono_VT, 'all'), sum(abs_dom_VT, 'all'));
fprintf(fid,'\nAverage Imp:%43.1f%28.1f%27.1f%26.1f', exact_star_average_imp, approx_star_average_imp, approx_zono_average_imp, abs_dom_average_imp);
fprintf(fid,'\nSafe Nets:%5d/%2d%18d/%2d%24d/%2d%22d/%2d%24d/%2d',Reluplex_safes, 45, exact_star_safes, 45, approx_star_safes, 45, approx_zono_safes, 45, abs_dom_safes, 45);
fclose(fid);
% 

% generate latex file for paper
% 
fid1 = fopen('P4_tab.tex', 'wt');
for i=1:5
    for j=1:9
        fprintf(fid1, '\n$N_{%d%d}$ & %s & %.2f & %s & %.2f & $%s{%.0f%s}$ & %s & %.2f & $%s{%.0f%s}$ & %s & %.2f & $%s{%.0f%s}$ & %s & %.2f & $%s{%.0f%s}$ %s', i, j, Reluplex_Res{i,j}, Reluplex_VT(i,j),exact_star_Res{i,j}, exact_star_VT(i,j),'\mathbf',Reluplex_VT(i,j)/exact_star_VT(i,j), '\times', approx_star_Res{i,j}, approx_star_VT(i,j),'\mathbf', Reluplex_VT(i,j)/approx_star_VT(i,j),'\times', approx_zono_Res{i,j}, approx_zono_VT(i,j),'\mathbf',Reluplex_VT(i,j)/approx_zono_VT(i,j),'\times', abs_dom_Res{i,j}, abs_dom_VT(i,j),'\mathbf',Reluplex_VT(i,j)/abs_dom_VT(i,j),'\times', '\\');
    end
end

fprintf(fid1,'\nTotal VT:$%s{%20.2f}$&$%s{%20.2f}$&$%s{%25.2f}$&$%s{%26.2f}$&$%s{%26.2f}$%s', '\mathbf', sum(Reluplex_VT, 'all'), '\mathbf', sum(exact_star_VT, 'all'), '\mathbf', sum(approx_star_VT, 'all'), '\mathbf', sum(approx_zono_VT, 'all'), '\mathbf', sum(abs_dom_VT, 'all'), '\\');
fprintf(fid1,'\nAverage Imp:$%s{%43.0f%s}$ & $%s{%28.0f%s}$ & $%s{%27.0f%s}$ & $%s{%26.0f%s}%s$', '\mathbf', exact_star_average_imp, '\times', '\mathbf', approx_star_average_imp, '\times', '\mathbf', approx_zono_average_imp, '\times', '\mathbf', abs_dom_average_imp, '\times', '\\');
fprintf(fid1,'\nSafe Nets:$%s{%5d/%2d}$ & $%s{%18d/%2d}$ & $%s{%24d/%2d}$ & $%s{%22d/%2d}$ & $%s{%24d/%2d}$%s','\mathbf',Reluplex_safes, 45, '\mathbf', exact_star_safes, 45, '\mathbf', approx_star_safes, 45, '\mathbf', approx_zono_safes, 45, '\mathbf', abs_dom_safes, 45, '\\');
fclose(fid);



