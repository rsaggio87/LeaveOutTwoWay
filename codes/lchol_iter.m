function Lchol=lchol_iter(xx)
    
%Credits to Chris Moser: https://github.com/rsaggio87/LeaveOutTwoWay/issues/17
    droptol=1e-2;
    diagcomp_linear_n = 1;
    diagcomp_linear_step = .1;
    diagcomp_linear_list = diagcomp_linear_step.*(1:diagcomp_linear_n);
    diagcomp_max = max(sum(abs(xx), 2)./diag(xx)) - 2; % value recommended by MATLAB to guarantee successful execution of -ichol()-, see https://www.mathworks.com/help/matlab/ref/ichol.html
    diagcomp_candidate_n = 20;
    diagcomp_candidate_base = 1.5;
    diagcomp_overshoot = 3;
    diagcomp_factor = 1./(diagcomp_candidate_base.^(diagcomp_candidate_n:-1:-diagcomp_overshoot));
    diagcomp_candidates = diagcomp_max.*diagcomp_factor;
    diagcomp_exponential_list = diagcomp_candidates(diagcomp_candidates > diagcomp_linear_step*diagcomp_linear_n);
    diagcomp_list = [diagcomp_linear_list diagcomp_exponential_list];
    Lchol = [];
    for diagcomp = diagcomp_list
        try
            Lchol=ichol(xx,struct('type','ict','droptol',droptol,'diagcomp',diagcomp));
            %fprintf('\n')
            %disp(['NOTE: function -ichol()- with diagcomp = ' num2str(diagcomp) ' succeeded!'])
            break % exit for loop after successful evaluation of -ichol()-
        catch
            disp(['USER WARNING: function -ichol()- with diagcomp = ' num2str(diagcomp) ' failed!'])
            if diagcomp == diagcomp_list(end)
                fprintf('\n')
                disp(['USER WARNING: function -ichol()- did not execute successfully for any value of diagcomp <= ' num2str(diagcomp)])
            end
        end
    end
end

