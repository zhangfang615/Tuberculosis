import pandas

if __name__ == '__main__':
    data = pandas.read_excel('E:/PYTHON_PROJECTS/TB_simulation/statistics.xlsx')
    output_file= 'E:/PYTHON_PROJECTS/TB_simulation/pearson_correlation.xlsx'
    data.corr().iloc[4:,0:4].to_excel(output_file)

    # ss = str(data.corr()).to_excel(output_file)
    # s1 = '\t'.join(filter(lambda x: x, ss.split(' ')))
    # output.write(s1)
    # output.close()