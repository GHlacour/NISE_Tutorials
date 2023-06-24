%% import_pdb.m
% * This function import .pdb files into the atom struct
%
%% Updated
% 20220811 by Long
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom('molecule.pdb','ABCDEF')
function atom = import_pdb(filename,chainlist)

% See http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
% 1 -  6         Record name   "ATOM  "
% 7 - 11         Integer       Serial       Atom  serial number.
% 13 - 16        Atom          Atom type    Atom name.   ->17 by MH
% 17             Character     AltLoc       Alternate location indicator.
% 18 - 20        Residue name  ResName      Residue name.
% 22             Character     ChainID      Chain identifier.
% 23 - 26        Integer       ResSeq       Residue sequence number.
% 27             AChar         Code         Code for insertion of residues.
% 31 - 38        Real(8.3)     X            Orthogonal coordinates for X in Angstroms.
% 39 - 46        Real(8.3)     Y            Orthogonal coordinates for Y in Angstroms.
% 47 - 54        Real(8.3)     Z            Orthogonal coordinates for Z in Angstroms.
% 55 - 60        Real(6.2)     Occupancy    Occupancy.
% 61 - 66        Real(6.2)     TempFactor   Temperature  factor.
% 73 - 76        LString(4)    Segment identifier, left-justified. % Not used
% 77 - 78        LString(2)    Element      Element symbol, right-justified.
% 79 - 80        LString(2)    Charge       Charge on the atom.
fid = fopen(filename,'r');
% fullText = fread(fid,'char=>char')';
% data = strread(fullText,'%s','delimiter','\n');% use textscan instead?
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); % New addition
data=data{1}; % New addition
fclose(fid);

IndexCRYS = strfind(data,'CRYS');
Index = find(not(cellfun('isempty',IndexCRYS)));

j = 0;atom=[];
for i = 1:length(data)
    line = data{i};
    if length(line)>=22 && contains(chainlist,line(22)) && strcmp(line(1:6),'HETATM')
        j = j + 1;
        atom(j).molid = str2double(line(23:26));
        atom(j).resname = strtrim(line(18:20));
        atom(j).chain = strtrim(line(22));
%         atom(j).type = {strtrim(line(13:16))}; % Changed to 17 for better
%         compatiblity
%         atom(j).fftype = {strtrim(line(13:16))}; % Changed to 17 for better
%         compatiblity
        atom(j).type = strtrim(line(13:17));
%         atom(j).fftype = {strtrim(line(13:17))};
        atom(j).index = str2double(line(7:11));
%         atom(j).neigh.type = {};
%         atom(j).neigh.index = [0;0;0;0;0;0];
%         atom(j).neigh.dist = [0;0;0;0;0;0];
%         atom(j).bond.type = [0;0;0;0;0;0];
%         atom(j).bond.index = [0;0;0;0;0;0];
%         atom(j).angle.type = [0;0;0;0;0;0];
%         atom(j).angle.index = [0;0;0;0;0;0];
        atom(j).x = str2double(line(31:38));
        atom(j).y = str2double(line(39:46));
        atom(j).z = str2double(line(47:54));
%         atom(j).vx = NaN;
%         atom(j).vy = NaN;
%         atom(j).vz = NaN;
        
%         occupancy(j,1)=str2double(line(55:60));
%         tempfactor(j,1)=str2double(line(61:66));
        
%         atom(j).occupancy=occupancy(j,1);
%         atom(j).B=tempfactor(j,1);
    end
end

nAtoms=size(atom,2);

fprintf('.pdb file imported from %s\n',filename)
