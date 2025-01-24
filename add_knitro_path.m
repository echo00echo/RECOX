if ispc
    addpath('C:\Program Files\Artelys\knitro-12.3.0-Win-64\knitromatlab');
elseif isunix
    addpath(strcat(getenv('KNITRO_HOME'),'/knitromatlab'));
elseif ismac
    error('add_knitro_path.m has not supported MacOS yet');
end
