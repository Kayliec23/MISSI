%%%
%%% createPBSfile.m
%%%
%%% Writes out a script that runs the MITgcm. Takes as arguments the name
%%% of the experiment, the number of nodes that will be used, and the wall
%%% time in hours (after which the computation will be cut off regardless 
%%% of whether it is complete).
%%%
function createPBSfile_Stampede(dirname,expname,nodes,walltime,account,cluster_path)
  %%% Limit walltime
  walltime = min(ceil(walltime),48);
  
  %%% Set variables correctly for PBS
  headertext = [...
  '#!/bin/bash \n' ...
  '#SBATCH -J ',expname,'                  # job name \n' ...
  '#SBATCH -o output_%%j.txt                # output and error file name (%%j expands to jobID)\n' ...
  '#SBATCH -n ',num2str(nodes),'           # total number of mpi tasks requested\n' ...
  '#SBATCH --nodes=',num2str(ceil(nodes/68)),'     # total number of nodes requested\n' ...
  '#SBATCH -p icx                       # queue (partition) -- normal, development, etc.\n' ...
  '#SBATCH -t ',num2str(walltime),':00:00  # run time (hh:mm:ss)\n' ...
  '#SBATCH --mail-user=kayliec23@g.ucla.edu\n' ...
  '#SBATCH --mail-type=begin               # email me when the job starts\n' ...
  '#SBATCH --mail-type=end                 # email me when the job finishes\n' ...  
  '\n' ...
  'cd ',cluster_path,'\n' ...
  'ibrun ./mitgcmuv']; ...
  %%% 'mpirun_rsh -np ',num2str(nodes),' a.out ./mitgcmuv'];
  %%% Open template script
  templatename = './DEFAULTS/results/run_mitgcm';
  tfid = fopen(templatename,'r');
  if (tfid == -1)
    error(['Could not read template SLURM run file: ',templatename]);
  end
  %%% Open output script and write header text
  wfid = fopen(fullfile(dirname,'run_mitgcm'),'w');
  if (wfid == -1)
    error('Could not open PBS script file');
  end
  fprintf(wfid,headertext); 
  %%% Copy over all template text
  count = 0;
  while true
    count = count+1;
    nextline = fgetl(tfid);    
    if ~ischar(nextline)
      break;
    end
    if (count <= 5) %%% skip placeholder comments at start of template
      continue;
    end
    fprintf(wfid,'%s\n',nextline);    
  end
  %%% Close files when we're done
  fclose(wfid);
  fclose(tfid);
end
