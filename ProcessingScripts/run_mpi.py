from mpi4py import MPI
import subprocess
import uuid
import os


# Common MPI communication, rank, size
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



protein = '3CLPro_protein'
smilefile = ''
pocket = '-10.520,-2.322,-20.631'
grid = '54,52,60'
residues = ''


# Which jobs have reported finishing at least once
reported = [0] * size

if rank == 0:

    big_smile_file = open('ena+db-small.can')
    smile_lines = big_smile_file.readlines()
    def get_next_large ():
        if len(smile_lines) > 0:
            line = smile_lines.pop().strip()
            return line
        return ''
    def get_next_small ():
        if len(smile_lines) > 0:
            line = smile_lines.pop(0).strip()
            return line
        return ''


    # Receive hello from worker and give instructions
    for i in range(1, size):
        datain = comm.recv(source=MPI.ANY_SOURCE)
        if i <= 1:
            dataout = {'line': get_next_large()}
        else:
            dataout = {'line': get_next_small()}
        comm.send(dataout, dest=datain['rank'])

    while True:
        datain = comm.recv(source=MPI.ANY_SOURCE)
        if datain['result'] < 0:
            break
        reported[datain['rank']] = 1
        if sum(reported) < size - 1:
            dataout = {'line': get_next_small()}
            comm.send(dataout, dest=datain['rank'])
        else:
            dataout = {'line': ''}
            comm.send(dataout, dest=datain['rank'])
            break

    for i in range(2, size):
        datain = comm.recv(source=MPI.ANY_SOURCE)
        dataout = {'line': ''}
        comm.send(dataout, dest=datain['rank'])

else:
    data = {'rank': rank}
    comm.send(data, dest=0)
    while True:
        data = comm.recv(source=0)
        if len(data['line']) == 0:
            break

        # Process line
        fname = '/tmp/'+str(uuid.uuid4())
        with open(fname, 'w') as fo:
            fo.write(data['line'] + '\n')
            print('myline', data['line'])

        print(['./smiles_dock.sh', protein, fname, pocket, grid])
        process = subprocess.run(['./smiles_dock.sh', protein, fname, pocket, grid])
        process = subprocess.run(['./summarize.sh', fname])

#        os.remove(fname)
        data = {'rank': rank, 'result': rank}
        comm.send(data, dest=0)




