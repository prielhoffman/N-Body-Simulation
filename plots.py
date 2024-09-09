import matplotlib.pyplot as plt

def scatter_data_from_files(file1, file2, file3):
    # Read data from the text files
    data1 = [] # start positions data
    with open(file1, 'r') as f1:
        for line in f1:
            x, y = map(float, line.strip().split())
            data1.append((x, y))

    data2 = [] # middle positions data
    with open(file2, 'r') as f2:
        for line in f2:
            x, y = map(float, line.strip().split())
            data2.append((x, y))

    data3 = [] # end positions data
    with open(file3, 'r') as f3:
        for line in f3:
            x, y = map(float, line.strip().split())
            data3.append((x, y))

    # scatter plot the data
    plt.scatter(*zip(*data1), label='start positions')
    plt.scatter(*zip(*data2), label='middle positions')
    plt.scatter(*zip(*data3), label='end positions')

    # add labels and legend
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()

    # show the plot
    plt.show()

# took the data from the output files
scatter_data_from_files('positions_start.txt', 'positions_middle.txt', 'positions_end.txt')
