function writeinput(niseStruct,fname)
% Write the structure niseStruct to the file fname
fid = fopen(fname,'w');
niseField = fieldnames(niseStruct);
niseValue = struct2cell(niseStruct);
for i = 1:length(niseField)
    fprintf(fid,'%s %s\n',niseField{i},num2str(niseValue{i}));
end
fclose(fid);
fprintf('%s generated\n',fname);
end