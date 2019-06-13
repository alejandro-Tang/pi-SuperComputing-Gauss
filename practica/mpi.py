import os

def main():
    clean_comand = 'make clean'
    compile_comand = 'make'
    os.system(clean_comand)
    os.system(compile_comand)
    for i in range(1, 4):
        os.system('mpirun -np ' + str(2^i) + ' gaussian matrix.2400.txt')

main()

