import os
import sys
import matplotlib.pyplot as plt

def main():
    clean_comand = 'make clean'
    compile_comand = 'make'
    os.system(clean_comand)
    os.system(compile_comand)
    test_time = []
    test_size = [8, 800, 1600, 2400, 3200]
    test_cores = []
    file_name = sys.argv[1]
    print("File name: ", file_name)

    for i in range(0, 4):
        out_time = os.popen('mpirun -np ' + str(2**i) + ' gaussian '+ file_name).read()
        test_cores.append(2**i)
        test_time.append(float(out_time[:-2]))
    print(test_time)
    plt.plot(test_cores,test_time)
    plt.show()
    pause = input()


main()

