from mpi4py import MPI
import subprocess
import uuid
import os
import time
import sys

# Common MPI communication, rank, size
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

protein = '3CLPro_protein'
smilefile = ''
pocket = '-10.520,-2.322,-20.631'
grid = '54,52,60'
residues = ''

smiles_file_name = sys.argv[1]

# Which jobs have reported finishing at least once
reported = [0] * size

if rank == 0:

    # Smile file passed in as argument
    big_smile_file = open(smiles_file_name)
    smile_lines = big_smile_file.readlines()

    # Methods to get lines from beginning and end
    def get_next_large ():
        if len(smile_lines) > 0:
            line = smile_lines.pop(0).strip()
            return line
        return ''
    def get_next_small ():
        if len(smile_lines) > 0:
            line = smile_lines.pop().strip()
            return line
        return ''

    # Record files with bad return codes
    fbad = open(smiles_file_name+'.bad', 'w')
    def record_bad (line):
        fbad.write(line+'\n')
        fbad.flush()
        return


    # Receive hello from worker and give instructions
    for i in range(1, size):
        datain = comm.recv(source=MPI.ANY_SOURCE)
        if i <= 1:
            dataout = {'line': get_next_large()}
        else:
            dataout = {'line': get_next_small()}
        comm.send(dataout, dest=datain['rank'])

    # delay in case some init bad
    time.sleep(5)

    # Go until there are no more
    while True:

        # Receive data from worker
        datain = comm.recv(source=MPI.ANY_SOURCE)
        if datain['result'] < 0:
            break
        if datain['result'] == 0:
            record_bad(datain['line'])

        # Record which workers have reported at least once
        reported[datain['rank']] = 1

        # Get the next smallest smile.  If next line is blank or number reported
        # is all-1 begin wrapping up
        nextline = get_next_small()
        if sum(reported) < size - 1 and len(nextline) > 0:
            dataout = {'line': nextline}
            comm.send(dataout, dest=datain['rank'])
        else:
            dataout = {'line': ''}
            comm.send(dataout, dest=datain['rank'])
            break

    # wait for workers to finish
    for i in range(2, size):
        datain = comm.recv(source=MPI.ANY_SOURCE)
        if datain['result'] == 0:
            record_bad(datain['line'])
        dataout = {'line': ''}
        comm.send(dataout, dest=datain['rank'])

else:
    # Initial data send from worker (just a hello)
    data = {'rank': rank, 'result': 0}
    comm.send(data, dest=0)

    # Run as long as we're getting good data
    while True:

        # Get data from master
        data = comm.recv(source=0)
        line = data['line']
        if len(line) == 0:
            print('done rank', rank)
            break

        # Process line by putting it in a temp file
        fname = '/tmp/'+str(uuid.uuid4())
        with open(fname, 'w') as fo:
            fo.write(line + '\n')

        # Process according to scripts
        print(['./smiles_dock.sh', protein, fname, pocket, grid])
        process1 = subprocess.run(['./smiles_dock.sh', protein, fname, pocket, grid])
        if process1.returncode == 0:
            process2 = subprocess.run(['./summarize.sh', fname])

        # Remove temp file
        os.remove(fname)

        # If there is an error detected return 0
        if process1.returncode != 0 or process2.returncode != 0:
            print('ERROR in:', line)
            data = {'rank': rank, 'line': line, 'result': 0}
        else:
            data = {'rank': rank, 'line': line, 'result': 1}
        comm.send(data, dest=0)




