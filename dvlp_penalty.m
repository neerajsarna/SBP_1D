function [penalty] = dvlp_penalty(An,B)
modAn = compute_Amod(An);
Xminus = compute_Xminus(An);
penalty = 0.5 * (An-modAn) * Xminus * inv(B*Xminus);
end