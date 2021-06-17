import pandas as pd
import numpy as np
from datetime import date
import os

class CRIC():
    def __init__(self):
        self.month_order = ['Januari', 'Februari', 'Maret', 'April', 'Mei', 'Juni', 'Juli', 'Agustus', 'September', 'October', 'November', 'Desember']
        self.stations = []
        self.wmoids = []
        self.latitudes = []
        self.longitudes = []
        self.altitudes = []
        self.inputfiles = []
        self.m_rrs = []
        self.m_tavgs = []
        self.n_stations = 0
    

    def read_excel(self, filepath):
        self.inputfile = pd.read_excel(filepath, 
                                  index_col='Tanggal', 
                                  skiprows=6, 
                                  na_values=['8888', ' ', '', '9999', '-9999'])
        self.inputfile.index = pd.to_datetime(self.inputfile.index, format='%d-%m-%Y')
        print(self.inputfile.head())
        meta = pd.read_excel(filepath, usecols='A:B', header=None)
        self.WMOID = meta[1][0]
        self.STATION = meta[1][1]
        self.LATITUDE = meta[1][2]
        self.LONGITUDE = meta[1][3]
        self.ALTITUDE = meta[1][4]
        print('WMOID :', self.WMOID)
        print('NAMA STASIUN :', self.STATION)
        print('lintang/buur :', self.LONGITUDE, '/', self.LATITUDE)
        print('Ketinggian :', self.ALTITUDE)
        print('Jumlah data kosong / jumlah data (Curah Hujan): ', self.inputfile['RR'].isna().sum(), '/', len(self.inputfile['RR']))
        print('Jumlah data kosong / jumlah data (Suhu udara rata-rata): ', self.inputfile['Tavg'].isna().sum(), '/', len(self.inputfile['Tavg']))
        # return self.inputfile
        self.daily2monnthly()
        self.add_to_stations()

    def add_to_stations(self):
        self.altitudes.append(self.ALTITUDE)
        self.ALTITUDE = None
        self.longitudes.append(self.LONGITUDE)
        self.LONGITUDE = None
        self.latitudes.append(self.LATITUDE)
        self.LATITUDE = None
        self.stations.append(self.STATION)
        self.STATION = None

        self.inputfiles.append(self.inputfile)
        self.inputfile = None

        self.dfRR = self.dfRR.rename('RR'+str(self.WMOID))
        self.m_rrs.append(self.dfRR)
        self.dfRR = None

        self.dfT = self.dfT.rename('Tavg'+str(self.WMOID))
        self.m_tavgs.append(self.dfT)
        self.dfT = None
        self.n_stations = self.n_stations + 1
        self.wmoids.append(self.WMOID)
        self.WMOID = None

        print('sudah ada '+str(self.n_stations)+ 'data stasiun')

    def daily2monnthly(self, S=[], agg_f=''):
        self.dfRR = self.inputfile['RR'].resample('M').sum()
        self.dfT = self.inputfile['Tavg'].resample('M').mean()

        # if agg_f == 'sum':
        #     df_tmp = S.resample('M').sum()
        # elif agg_f == 'mean':
        #     df_tmp = S.resample('M').mean()
        # return df_tmp

    def CDD(self,S):
    #  print('Shape S: ', S.shape)
        
        ind_CDD=[]
        S_no_nan = S[~np.isnan(S)]
        N = len(S)
        N2 = len(S_no_nan)
        if ((N2/N) < 0.3): 
            ind_CDD = np.empty(1)
            ind_CDD = np.nan
        else:
            temp = 0
            ind_CDD = 0 
            j =0
            while (j < N2):
                while (j < N2 ) and (S_no_nan[j] < 1.0 ):
                    j += 1
                    temp +=1
                if ind_CDD < temp:
                    ind_CDD = temp
                temp = 0
                j += 1 
        return ind_CDD
 
    def CWD(self,S):
        ind_CWD=[]
        S_no_nan = S[~np.isnan(S)]
        N = len(S)
        N2 = len(S_no_nan)
        if ((N2/N) < 0.3): 
            ind_CWD = np.empty(1)
            ind_CWD = np.nan
        else:
            temp = 0
            ind_CWD = 0 
            j =0
            while (j < N2):
                while (j < N2 ) and (S_no_nan[j] > 1.0 ):
                    j += 1
                    temp +=1
                if ind_CWD < temp:
                    ind_CWD = temp
                temp = 0
                j += 1 
        return ind_CWD
 
 
    def rx1day(self,S):
        return S.max()
    
    def rx3day(self, S):
        ind_R3d=[]
        S_no_nan = S[~np.isnan(S)]
        N = len(S)
        N2 = len(S_no_nan)
        if ((N2/N) < 0.3): 
            ind_R3d = np.empty(1)
            ind_R3d = np.nan
        else:
            temp = 0
            ind_R3d = 0 
            for i in range(0,N-2):
                if (~np.isnan(S[i])) and  (~np.isnan(S[i+1]))  and  (~np.isnan(S[i+2])):
                    temp = S[i] + S[i+1] + S[i+2]
                if ind_R3d < temp:
                    ind_R3d = temp
        return ind_R3d
  
    def rx5day(self, S):
        ind_R3d=[]
        S_no_nan = S[~np.isnan(S)]
        N = len(S)
        N2 = len(S_no_nan)
        if ((N2/N) < 0.3): 
            ind_R5d = np.empty(1)
            ind_R5d = np.nan
        else:
            temp = 0
            ind_R5d = 0 
            for i in range(0,N-4):
                if (~np.isnan(S[i])) and  (~np.isnan(S[i+1]))  and  (~np.isnan(S[i+2]))  and  (~np.isnan(S[i+3])) and  (~np.isnan(S[i+4])):
                    temp = S[i] + S[i+1] + S[i+2] + S[i+3] + S[i+4]
                if ind_R5d < temp:
                    ind_R5d = temp
        return ind_R5d
 
    def r20mm(self, S):
        return np.count_nonzero(S >= 20)
    
    def Prec95p(self, S):
        return S.quantile(q=0.95)

    def calculateIndices(self, S):
        self.data_cdd = self.inputfile['RR'].resample('Y').agg(C.CDD)
        self.data_cwd = self.inputfile['RR'].resample('Y').agg(C.CWD)
        self.data_rx1day = self.inputfile['RR'].resample('Y').agg(C.rx1day)
        self.data_rx3day = self.inputfile['RR'].resample('Y').agg(C.rx3day)
        self.data_rx5day = self.inputfile['RR'].resample('Y').agg(C.rx5day)
        self.data_r20mm = self.inputfile['RR'].resample('Y').agg(C.r20mm)
        self.data_Prec95p = self.inputfile['RR'].resample('Y').agg(C.Prec95p)

    # def generate_sibias_input(self, S, outfile):
    #     print('SiBiaS input untuk '+str(self.n_stations) + ' akan dibuat')

    #     ########## generate global timerange


    #     first_col = ['Lat', 'Lon'] + ["01/%02d/%04d"%(i.month, i.year) for i in S.index.to_pydatetime()]
    #     # first_col = ['Lat', 'Lon'] + [date(i.year, i.month, 1) for i in S.index.to_pydatetime()]
    #     outdict = {}
    #     outdict['Time'] = first_col
    #     second_col = [str(self.LATITUDE), str(self.LONGITUDE), ] + S.to_list()
    #     outdict['ID'+str(self.WMOID)] = second_col

    #     dfo = pd.DataFrame(outdict)
    #     # outfile='prec_'+str(C.WMOID)+'.xlsx'
    #     pd.io.formats.excel.header_style = None
    #     dfo.to_excel(outfile, index=False)
    
    def generate_sibias_input(self,outname='output.xlsx', var='RR'):
        mins = []
        maxs = []

        if self.n_stations < 4:
            raise 'Jumlah stasiun kurang, Tambah stasiun hingga minimal 4 stasiun'

        if var == 'RR':
            S = self.m_rrs
        elif var == 'Tavg':
            S = self.m_tavgs
        else:
            raise 'Variabel tidak dikenal'
        
        for i in range(len(S)):
            # dates.append(C.m_rrs[i].index)
            mins.append(S[i].index.min())
            maxs.append(S[i].index.max())
        index_min = np.array(mins).min()
        index_max = np.array(maxs).max()

        # var='RR'
        dfout = pd.DataFrame(index=pd.date_range(start=index_min, end=index_max, freq='M'))
        # dfout.append()
        dfout = pd.concat(S, axis=1)# dfout[C.m_rrs[0]]
        dfout = dfout.fillna(-9999)
        # dfout.index = dfout.index.map(lambda x: date(x.year, x.month, 1))
        # dfout.to_excel('4stationstest.xlsx')

        first_col = ['Lat', 'Lon'] + ["01/%02d/%04d"%(i.month, i.year) for i in dfout.index.to_pydatetime()]
        outdf = {}
        outdf['Time'] = first_col
        for i in range(self.n_stations):
            col = [self.latitudes[i], self.longitudes[i]] + dfout[var+str(self.wmoids[i])].to_list()
            outdf[var+str(self.wmoids[i])] = col
        
        pd.DataFrame(outdf).to_excel(outname, index=False)
        print('File '+outname+' telah berhasil dibuat')

    def read_excel_bmkg(self, filepath='/content/drive/MyDrive/Stasiun Meterologi Japura'):
        years = sorted([int(i) for i in os.listdir(filepath)])
        all_data = pd.DataFrame()
        for y in years:
            m_files = os.listdir(os.path.join(filepath, str(y)))
            for m in self.month_order:
                print('y , m', y, m)
                if os.path.exists(os.path.join(filepath, str(y), m+'_'+str(y)+'_laporan_iklim_harian.xlsx')):
                    print('turee')
                    df_tmp = pd.read_excel(os.path.join(filepath, str(y), m+'_'+str(y)+'_laporan_iklim_harian.xlsx'), 
                                  index_col='Tanggal',
                                  usecols= 'A,D,F',
                                  nrows=31, 
                                  skiprows=8, 
                                  na_values=['8888', ' ', '', '9999', '-9999'])
                    df_tmp = df_tmp.dropna(how='all')
                    df_tmp.index = pd.to_datetime(df_tmp.index, format='%d-%m-%Y')
                else:
                
                    df_tmp = pd.DataFrame(
                        {'Tanggal' : pd.date_range(start='%02d/%4d'%(month_order.index(m)+1, y), 
                                                end = ('%02d/%4d'%(month_order.index(m)+2, y)) if (month_order.index(m)+2 != 13) else ('%02d/%4d'%(1, y+1)))[:-1],
                        'Tavg' : np.nan,
                        'RR': np.nan
                        },
                    )
                    df_tmp = df_tmp.set_index('Tanggal')
                all_data = all_data.append(df_tmp)
        return all_data