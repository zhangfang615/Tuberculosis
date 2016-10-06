import pandas

if __name__ == '__main__':
    for j in range (1,8):
        output_file = 'E:/PYTHON_PROJECTS/TB_simulation/statistics_'+str(j)+'/8-12.xlsx'
        data = pandas.read_excel('E:/PYTHON_PROJECTS/TB_simulation/statistics_'+str(j)+'/simulation_statistics_10.xlsx')
        data = data[(data['resistant_ratio'] >= 8) & (data['resistant_ratio'] <= 12)]
        for i in range(2, 6):
            D = pandas.read_excel(
                'E:/PYTHON_PROJECTS/TB_simulation/statistics_'+str(j)+'/simulation_statistics_' + str(i) + '0.xlsx')
            filtered = D[(D['resistant_ratio'] >= 8) & (D['resistant_ratio'] <= 12)]
            data = data.append(filtered)
        data = data.reset_index(drop=True)
        data.loc[:, ['parameters', 'resistant_ratio']].to_excel(output_file)
    # data.to_excel(output_file)