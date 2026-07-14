% This file imports all coordinates from the created text file to COORD
clear all
clc
name = '2024-09-13_R_40x1xlsm5_plast_viso0_49_UniAx_0_COORD';

opts = detectImportOptions(name);
opts = setvartype(opts, 'double');
% Deaktiviere das Erkennen von Spaltennamen aus der Datei
% Spezifiziere, dass keine Variablennamen aus der Datei gelesen werden
opts = setvaropts(opts, "Part_Instance", "TrimNonNumeric", true);
opts = setvaropts(opts, "Part_Instance", "ThousandsSeparator", ",");

COORD = readtable([name,'.txt'], opts);
COORD(1:3, :) = [];
save([name(1,12:end),'.mat'], 'COORD');