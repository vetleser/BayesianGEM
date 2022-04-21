from mpi4py import MPI

comm=MPI.COMM_WORLD
myrank=comm.Get_rank()
ranksize=comm.Get_size()
print("hello world from rank {}, of {} ranks".format(myrank,ranksize))