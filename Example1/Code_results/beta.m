


A = beta_save;
% Save as an Excel file
%for kk = 1:size(A,3)
for kk = 1:100
xlswrite('n1000p4500_beta_p_ar1.xlsx',A(:,:,kk),kk);
end


B =  double(samp_save);
% Save as an Excel file
%for kk = 1:size(B,3)
for kk = 1:100
  xlswrite('n1000p4500_sample_p_ar1.xlsx',B(:,:,kk),kk);
end

C =  double(beta_supp_est_save);
% Save as an Excel file
%for kk = 1:size(C,3)
for kk = 1:100
  xlswrite('n1000p4500_screening_sample_p.xlsx',C(:,:,kk),kk);
end

D =  double(betO_save);
% Save as an Excel file
writematrix(D, "Oracle_1.xlsx")





