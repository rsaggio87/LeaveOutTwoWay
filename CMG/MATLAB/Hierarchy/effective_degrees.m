% function deg = effective_degrees(A)
%
% deg is the effective degree sequence of the adjacency matrix A: 
%  max_incident_weight/total_incident_weight

function deg = effective_degrees(A)


deg = full(max(A)./sum(A)); 


