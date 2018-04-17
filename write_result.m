function write_result(result,filename)
dlmwrite(filename,result(1).X','delimiter','\t','precision',10);

for j = 1 : length(result)
    dlmwrite(filename,result(j).sol','-append','delimiter','\t','precision',10);
end

end