%%%
%%% write_data_kpp
%%%
%%% Writes the 'data.kpp' input file from 
%%% the cell array KPP_PARM of parmlist objects
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_ggl90 (dirname,GGL90_PARM,listterm,realfmt)

  %%% Open the 'data.kpp' file for writing
  fname = 'data.ggl90';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# ===================\r\n' ...
    '# | GGL90 parameters  |\r\n' ...
    '# ===================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);
  
    %%% Parameter section titles
  titles = {...
    'Mixing length limit', ...
    'Minimum surface boundary condition', ...
    };
  
    %%% For each parameter section 
  for nparam=1:1:length(GGL90_PARM)
    
    %%% Write section header
    fprintf(fid,['# ',titles{nparam},'\r\n']);
    fprintf(fid,[' &GGL90_PARM0',num2str(nparam),'\r\n']);   
        
    %%% Write each parameter out to the 'data' file
    nextparmlist = GGL90_PARM{nparam};
    for n=1:1:nextparmlist.getLength()      
      writeParam(fid,nextparmlist.getParm(n),realfmt);
    end    
    
    %%% Write section footer
    fprintf(fid,[' ',listterm,'\r\n']);
    fprintf(fid,'\r\n');
    
  end
  
  %%% Close the file when we're finished
  fclose(fid);

end