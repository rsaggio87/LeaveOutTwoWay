% function [cI, nc, D, C] = steiner_group(A,B)
%
% Input A: adjacency matrix A
%          [optional] Matrix B with at most one not zero entry per column 
%        
%
% Output: Vertex disjoint clusters given in cI
%         Number of clusters nc
%         Vector D containing vertex volumes  
%         Matrix C with at most one not zero entry per column - The
%          connected components of C are the disjoint clusters in cI
%         Input B is used as a first approximation to matrix C


%/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
%/* 
%/* The CMG solver is distributed under the terms of the GNU General Public   */
%/* Lincense Version 3.0 of the Free Software Foundation.                    */
%/* CMG is also available under other licenses; contact authors for details. */


function [cIout, nc, D, C] = steiner_group(A,B)

if (nargin == 1)
    [C1,D,M] = graphprofile(A); 
    C = splitforest(C1);
    eff_deg = M./D;
    nodes_to_check = find(eff_deg<(1/8)); %ad-hoc set to 1/16
    if (isempty(nodes_to_check)==0)
        [C1,D] = update_groups(A,C);
        C=C1;
    end
    [cIout nc] = forest_components(C); 
else
    n = length(A);
    [C,D] = update_groups(A,B);
    [cIout nc] = forest_components(C);

    if (nc/n)>0.7
        [cIout, D,C] = steiner_group(A);
    end
end



