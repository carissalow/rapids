import numpy as np
import math
import pandas as pd
from scipy.stats import *
from tqdm import tqdm
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist, squareform
import datetime
import glob
from dateutil import tz


np.set_printoptions(suppress=True, formatter={'float_kind':'{:f}'.format})
np.random.seed(22)

def run_test(output_file):
    print("Output DATA: ")
    print(output_file)

def run_barnett_features_for_rapids(input_dataframe, accuracy_limit=51.0, timezone="", wtype="GLR", spread_pars=[10,1], center_rad=200, interval=10, n_reps=1, min_pause_dur=300, min_pause_dist=60, r=None, w=None, tint_m=None, tint_k=None):
    data_frame = input_dataframe
    lonlat, r, w = preprocessing(data_frame, interval=interval, acc_threshold=accuracy_limit, r=r, w=w, tint_m=tint_m, tint_k=tint_k)
    mobmatmiss = convert_to_flights_pauses(lonlat, r, w)
    mobmat = guess_pause(mobmatmiss, min_pause_dur, min_pause_dist)
    obj = initialize_params(mobmat)
    data_for_pandas = []
    if obj:
        lsmf = []
        lssigloc = []
        for i in range(n_reps):
            print("Sim #: ", i+1)
            out3 = simulate_mobility_gaps(mobmat, obj, wtype, spread_pars)
            out3 = np.array(out3)
            IDundef=np.where(out3[:,0]==3)[0]

            if len(IDundef) > 0:            
                out3 = np.delete(out3, IDundef,axis=0)
            obj3 = initialize_params(out3)

            output_features, slout, row_names = get_mobility_features(out3, obj3, mobmatmiss, timezone, center_rad, interval)
            lsmf.append(output_features)
            lssigloc.append(slout)

        #get average
        np_lsmf = np.array((lsmf))
        avg_output_features = [] #convert to pandas dataframe and then return this
        for j in range(len(row_names)):
            temp_avg = []
            for i in range(n_reps):
                temp_avg.append(np_lsmf[i][j])
            avg_output_features.append(np.mean(temp_avg, axis=0))

        
        for idx, row in enumerate(avg_output_features):
            date = row_names[idx]
            date_str = str(date[0]) + "-" + str(date[1]) + "-" + str(date[2])
            features = [date_str]
            for feature in row:
                features.append(feature)
            data_for_pandas.append(features)

    column_names = ["local_date", "hometime", "disttravelled", "rog", "maxdiam", "maxhomedist", "siglocsvisited", "avgflightlen", "stdflightlen", "avgflightdur", "stdflightdur", "probpause", "siglocentropy", "minsmissing", "circdnrtn", "wkenddayrtn"]
    df = pd.DataFrame(data_for_pandas, columns = column_names)
    df.set_index('local_date', inplace=True)
    return df

def run_barnett_features(input_dir, output_file, wtype="GLR", spread_pars=[10,1], timezone="", center_rad=200, interval=10, acc_threshold=51.0, n_reps=1, min_pause_dur=300, min_pause_dist=60, r=None, w=None, tint_m=None, tint_k=None):
    data_frame = load_beiwe(input_dir)
    lonlat, r, w = preprocessing(data_frame, interval=interval, acc_threshold=acc_threshold, r=r, w=w, tint_m=tint_m, tint_k=tint_k)
    mobmatmiss = convert_to_flights_pauses(lonlat, r, w)
    mobmat = guess_pause(mobmatmiss, min_pause_dur, min_pause_dist)
    obj = initialize_params(mobmat) 
    lsmf = []
    lssigloc = []
    for i in range(n_reps):
        print("Sim #: ", i+1)
        out3 = simulate_mobility_gaps(mobmat, obj, wtype, spread_pars)
        out3 = np.array(out3)
        IDundef=np.where(out3[:,0]==3)[0]

        if len(IDundef) > 0:            
            out3 = np.delete(out3, IDundef,axis=0)
        obj3 = initialize_params(out3)

        output_features, slout, row_names = get_mobility_features(out3, obj3, mobmatmiss, timezone, center_rad, interval)
        lsmf.append(output_features)
        lssigloc.append(slout)

    #get average
    np_lsmf = np.array((lsmf))
    avg_output_features = []
    for j in range(len(row_names)):
        temp_avg = []
        for i in range(n_reps):
            temp_avg.append(np_lsmf[i][j])
        avg_output_features.append(np.mean(temp_avg, axis=0))

    print("Writing features to ", output_file)
    with open(output_file, "w") as writer:
        writer.write("Date\tHometime\tDistTravelled\tRoG\tMaxDiam\tMaxHomeDist\tSigLogsVisited\tAvgFlightLen\tStdFlightLen\tAvgFlightDur\tStdFlightDur\tProbPause\tSigLocEntropy\tMinsMissing\tCircdnRtn\tWkEndDayRtn\n")
        for idx, row in enumerate(avg_output_features):
            date = row_names[idx]
            date_str = str(date[0]) + "-" + str(date[1]) + "-" + str(date[2])
            feature_str = date_str
            for feature in row:
                feature_str += "\t"+str(feature)
            writer.write(feature_str + "\n")

    print("Finish writing to ", output_file)


def load_beiwe(input_dir, debug_mode=True):
    # load all files ending with .csv in input_dir    
    print("Reading csv files from {}".format(input_dir))
    all_files = glob.glob(input_dir + "/*.csv")
    the_list = []

    counter = 1 
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        the_list.append(df)
        counter += 1

    new_frame = pd.concat(the_list, axis=0, ignore_index=True)
    return new_frame

#Convert GPS raw data to mobile matrix (GPS2MobMat.R)
def preprocessing(df, interval=10.0, acc_threshold=51.0, r=None,w=None,tint_m=None,tint_k=None):
    
    df.columns = ['timestamp','latitude','longitude','altitude','accuracy']
    
    # ===============================
    # 1. order by the timestamp (ascending)
    # ===============================
    df = df.sort_values(by=['timestamp'])
    # print(df)
    
    # ===============================
    # 2. select row if accuracy < threshold #pass
    # ===============================
    df_acc = df[df['accuracy'] < acc_threshold]
    df = df_acc

    # ===============================
    #3. collapse data to interval #pass
    # ===============================
    if(tint_k != None and tint_m != None):
        t0 = np.ceil(df['timestamp'].values[0]/1000)
        df = df[df['timestamp']/1000 - t0 % (tint_k + tint_m) < tint_k]

    if not r:
        r = np.sqrt(interval)

    if not w:        
        w = np.mean(df['accuracy'].values) + interval

    start_time = (df['timestamp'].iloc[0]/1000)
    end_time = (df['timestamp'].iloc[-1]/1000)

    avg_matrix_lst = []
    
    IDam = 0
    count = -1
    next_line = [1, start_time + interval/2.0, df['latitude'].iloc[0], df['longitude'].iloc[0], 0.0, 0.0]    
    interval_counter = 1
    print("Collapse data within", interval, "second intervals...\n")

    counter = 0
    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        if counter > 0:
            timestamp = row['timestamp']
            latitude = row['latitude']
            longitude = row['longitude']

            if timestamp/1000.00 < start_time + interval:
                next_line[2] += latitude
                next_line[3] += longitude
                interval_counter += 1
                
            else:
                next_line[2] = next_line[2]/interval_counter
                
                next_line[3] = next_line[3]/interval_counter

                temp = next_line                
                avg_matrix_lst.append(temp)

                num_miss = np.floor((timestamp/1000 - (start_time + interval))/interval)

                if num_miss > 0:
                    has_miss = [4, start_time + interval/2.0, start_time + interval * (num_miss + 1) +interval/2.0, float("NaN"), float("NaN"), float("NaN")]
                    avg_matrix_lst.append(has_miss)

                                
                start_time = start_time + interval * (num_miss+1)
                next_line = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                next_line[0] = 1
                next_line[1] = start_time + interval/2.0
                next_line[2] = latitude
                next_line[3] = longitude
                interval_counter = 1
        counter += 1

    avg_matrix_df = pd.DataFrame(data=avg_matrix_lst)
    
    # ================================================================
    # avg_matrix_df is a matrix of average scores per interval
    # avg_matrix_df[0] = sign of the case, 1 = normal, 4 = has missing value
    # avg_matrix_df[1] = time
    # avg_matrix_df[2] = latitude
    # avg_matrix_df[3] = longitude
    # avg_matrix_df[4] = x
    # avg_matrix_df[5] = y
    # ================================================================
    
    # ===============================
    #4. convert from Lat/Lon to X/Y
    # ===============================        
    print("Convert latitude/longitude to X/Y\n")
    new_dataframe = lat_long_to_xy(avg_matrix_df)
    return new_dataframe, r, w

def convert_to_flights_pauses(new_dataframe, r, w, output_file=None):
    print("\nConvert from X/Y to flights/pauses...\n")
    
    mobmatmiss = extract_flights_from_dataframe(new_dataframe, r, w)
    

    if output_file: #output_file = mobmatmiss
        print("Writing {} rows to mobmatmiss to tsv...".format(mobmatmiss.shape[0]))
        mobmatmiss.to_csv(output_file, sep="\t")
        print("Finish writing mobmatmiss to tsv...")
    return mobmatmiss

def extract_flights_from_dataframe(dataframe, r, w):
    np_arr = dataframe.to_numpy()
    cur_idx = 0
    result_matrix = []
    temp_counter = 0
    prev_rows_counter = 0
    for i in tqdm(range(len(np_arr)), total=len(np_arr)):
        if np_arr[i][0] == 4:                        
            row_with_missing_data = []
            prev_rows = extract_flights(np_arr[cur_idx:i], r, w)            
            last_timestamp = np_arr[i][1]

            row_with_missing_data.append([4, float("NaN"), float("NaN"), last_timestamp, float("NaN"), float("NaN"), np_arr[i][2]])    
            result_matrix += prev_rows            
            result_matrix += row_with_missing_data
            prev_rows_counter += len(prev_rows)

            temp_counter += 1
            cur_idx = i+1

    if cur_idx < len(np_arr):
        rows = extract_flights(np_arr[cur_idx:len(np_arr)], r, w)
        last_timestamp = rows[len(rows)-1][6]
        result_matrix += rows
        
        if last_timestamp < np_arr[len(np_arr)-1][1]:
            if result_matrix[len(result_matrix)-1][0] != 2:
                pause = [2, rows[len(rows)-1][4], rows[len(rows)-1][5], last_timestamp, float("NaN"), float("NaN"), np_arr[len(np_arr)-1][1]]
                result_matrix.append(pause)
            else:
                result_matrix[len(result_matrix)-1][6] = np_arr[len(np_arr)-1][1]

    result_df = pd.DataFrame(data=result_matrix, index=None, columns=["code", "lon1", "lat1", "t1", "lon2", "lat2", "t2"])
    return result_df


# input1 = x_list --> np array [x0, x1, x2, ...xn]
# input2 = y_list
def is_flight(x_listt, y_listt, r, w):
    
    x_list = np.array(x_listt)
    y_list = np.array(y_listt)
    num_of_rows = len(x_list)
    prev_lat_lon = np.array((x_list[0], y_list[0])) #x0, y0
    next_lat_lon = np.array((x_list[num_of_rows-1], y_list[num_of_rows-1]))
    # euclidean_dist = euclidean distance of [x1,y1] and [x2, y2]
    single_euclidean_dist = np.linalg.norm(prev_lat_lon - next_lat_lon)
    if single_euclidean_dist < r: #DONE
        return False

    next_x = x_list[1:num_of_rows]
    prev_x = x_list[0:num_of_rows-1]
    subs_x = next_x - prev_x
    pow_x = np.power(subs_x, 2)

    next_y = y_list[1:num_of_rows]
    prev_y = y_list[0:num_of_rows-1]
    subs_y = next_y - prev_y
    pow_y = np.power(subs_y, 2)

    sum_sqrt_x_y = np.sqrt(pow_x + pow_y)

    if np.min(sum_sqrt_x_y) < r:
        return False 
    
    if num_of_rows == 2:
        return True

    if x_list[0] == x_list[num_of_rows-1]:
        if np.max(np.abs(x_list[1: num_of_rows])) > w: 
            return False
        else:
            return True

    if x_list[0] > x_list[num_of_rows-1]:
        x_list = x_list[::-1]
        y_list = y_list[::-1]

    x_new = x_list - x_list[0]
    y_new = y_list - y_list[0]
    theta = -1*np.arctan(y_new[num_of_rows-1]/x_new[num_of_rows-1])
    divvv = y_new[num_of_rows-1]/x_new[num_of_rows-1]
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    a_matrix =  np.array(([[costheta, -1*sintheta], [sintheta, costheta]]))
    temp_matrix = np.array((x_new[1:num_of_rows-1], y_new[1:num_of_rows-1]))
    rotation_pts = np.matmul(a_matrix, temp_matrix)

    
    if np.max(np.abs(rotation_pts[1])) > w:
        return False
    else:
        return True

def lat_long_to_xy(df,R=6371000.0): #1000.000
    df1 = df[df[0] == 1.0].copy()
    
    lat_v,lon_v = df1[2].values, df1[3].values
    th0 = np.min(lon_v)    
    th1 = np.max(lon_v)
    ph0 = np.min(lat_v)
    ph1 = np.max(lat_v)

    pi = math.pi
    d1 = 2*pi*R*((ph1-ph0)*2*pi/360)/(2*pi)
    d1 = 2*pi*R*((ph1-ph0)*2*pi/360)/(2*pi)
    d2 = 2*pi*(R*np.sin(pi/2-ph1*2*pi/360))*((th1-th0)*2*pi/360)/(2*pi)
    d3 = 2*pi*(R*np.sin(pi/2-ph0*2*pi/360))*((th1-th0)*2*pi/360)/(2*pi)
    
    for idx, row in df.iterrows():
        case = row[0]
        latitude = row[2]
        longitude = row[3]
        if case == 1:
            w1=(latitude-ph0)/(ph1-ph0)
            w2=(longitude-th0)/(th1-th0)
            # df[4][idx]= w1 * np.abs(d3-d2)/2+w2*(d3*(1-w1)+d2*w1)
            df.loc[idx, 4] = w1 * np.abs(d3-d2)/2+w2*(d3*(1-w1)+d2*w1)
            # df[5][idx]= w1 * d1 * np.sin(np.arccos(np.abs((d3-d2)/(2*d1))))    
            df.loc[idx, 5] = w1 * d1 * np.sin(np.arccos(np.abs((d3-d2)/(2*d1))))    
    
    return df


def extract_flights(matrix, r, w):
    np_matrix = np.array(matrix)
    if len(np_matrix) == 1:
        row = [[3, np_matrix[0][4], np_matrix[0][5], np_matrix[0][1], float("NaN"), float("NaN"), float("NaN")]]
        return row
    
    cur_idx = 0
    output = []
    timestamp = np_matrix[cur_idx][1] 
    x = np_matrix[cur_idx][4] 
    y = np_matrix[cur_idx][5] 

    distance_line = [x, y, timestamp, float("NaN"), float("NaN"), float("NaN")]
    flight_stat = False
    last_pause_time = 0
    last_pause_idx = -1
    while True:
        
        next_idx = cur_idx + 1
        if next_idx == len(np_matrix)-1:
            distance_line[3] = np_matrix[next_idx][4]
            distance_line[4] = np_matrix[next_idx][5]
            distance_line[5] = np_matrix[next_idx][1]
            output.append(distance_line)
            break
        
        while True:            
            part_of_matrix = np_matrix[cur_idx:next_idx+1]    
            x_s = part_of_matrix[:,4]
            y_s = part_of_matrix[:,5]
            
            flight_stat = is_flight(x_s, y_s, r, w)

            if not flight_stat:
                last_pause_time = part_of_matrix[len(part_of_matrix)-1][1]            
                last_pause_idx = next_idx
                break
            
            next_idx += 1

            if next_idx >= len(np_matrix):
                break


        if next_idx == cur_idx + 1 and next_idx < len(np_matrix):            
            distance_line[2] = np_matrix[next_idx][1]            
            np_matrix = np.delete(np_matrix, (next_idx), axis=0)
        else:
            distance_line[3] = np_matrix[next_idx-1][4]
            distance_line[4] = np_matrix[next_idx-1][5]
            distance_line[5] = np_matrix[next_idx-1][1]
            
            output.append(distance_line)            
            cur_idx = next_idx - 1

            distance_line = [np_matrix[cur_idx][4], np_matrix[cur_idx][5], np_matrix[cur_idx][1], float("NaN"), float("NaN"), float("NaN")]

        if next_idx >= len(np_matrix):
            break

    if len(output) == 0:
        last_timestamp_in_matrix = np_matrix[len(np_matrix)-1][1]
        row = [2, x, y, np_matrix[0][1], float("NaN"), float("NaN"), last_timestamp_in_matrix]
        return [row]

    result = []
    if timestamp < output[0][2]:        
        row = [2, x, y, timestamp, float("NaN"), float("NaN"), output[0][2]]
        result.append(row)
        if len(output) == 1:
            result[0][6] = output[0][5]
            return result

    if len(output) == 1:
        result = [[1] + output[0]]
        return result

    for idx, row in enumerate(output):
        if result != []:
            cur_time = row[2]            
            prev_row = result[len(result)-1]
            if prev_row[0] == 2: #is a pause
                if cur_time == prev_row[6]:
                    row[0] = prev_row[1]
                    row[1] = prev_row[2]
                    result.append([1] + row)
                else:
                    print("error sanity check") #sanity check
            else:
                if prev_row[6] != cur_time:
                    #add pause
                    prev_x = prev_row[4]
                    prev_y = prev_row[5]
                    prev_t = prev_row[6]
                    pause_row = [2, prev_x, prev_y, prev_t, float("NaN"), float("NaN"), cur_time]
                    result.append(pause_row)
                    flight_row = [1, prev_x, prev_y, cur_time, row[3], row[4], row[5]]
                    result.append(flight_row)
                else:
                    result.append([1] + row)
        else:            
            result.append([1] + row)

    last_row = result[len(result)-1]
    last_timestamp = last_row[6]
    result = np.array(result)
    ID_flight = np.where(result[:,0] == 1)[0]
    ID_pause = np.where(((result[ID_flight,1] - result[ID_flight,4]) ** 2 + (result[ID_flight,2]-result[ID_flight,5])**2) == 0)[0] 

    if len(ID_pause) > 0:
        for i in len(ID_pause):
            result[ID_flight[ID_pause], 0] = 2
            result[ID_flight[ID_pause], 4] = float("nan")
            result[ID_flight[ID_pause], 5] = float("nan")

    result = list(result)
    
    if last_timestamp < last_pause_time:
        status = last_row[0]
        if status == 1:
            x_0 = last_row[4]
            y_0 = last_row[5]
        else:
            x_0 = last_row[1]
            y_0 = last_row[2]
        x_1 = np_matrix[last_pause_idx][4]
        y_1 = np_matrix[last_pause_idx][5]
        flight_stat = is_flight(np.array([x_0, x_1]), np.array([y_0, y_1]), r, w)
        if flight_stat:                
            end_row = [1, x_0, y_0, last_timestamp, x_1, y_1, last_pause_time]        
        else:
            
            end_row = [2, x_0, y_0, last_timestamp, float("NaN"), float("NaN"), last_pause_time]
        result.append(end_row)
    
    last_row_in_matrix = np_matrix[len(np_matrix)-1]
    last_timestamp_in_matrix = last_row_in_matrix[1]

    last_row = result[len(result)-1]
    last_timestamp = last_row[6]
    
    if last_timestamp < last_timestamp_in_matrix:
        status = last_row[0]
        if status == 1:
            x_0 = last_row[4]
            y_0 = last_row[5]
        else:
            x_0 = last_row[1]
            y_0 = last_row[2]
        x_1 = last_row_in_matrix[4]
        y_1 = last_row_in_matrix[5]
        flight_stat = is_flight(np.array([x_0, x_1]), np.array([y_0, y_1]), r, w)

        if flight_stat:
            end_row = [1, x_0, y_0, last_timestamp, x_1, y_1, last_timestamp_in_matrix]
            result.append(end_row)
            
        else:
            if result[len(result)-1][0] == 2:
                result[len(result)-1][6] = last_timestamp_in_matrix
            else:
                end_row = [2, x_0, y_0, last_timestamp, float("NaN"), float("NaN"), last_timestamp_in_matrix]
                result.append(end_row)
    return result

def guess_pause(matrix, min_duration=300, pause_rad=75):
    print("Inferring pauses...")
    np_matrix = np.array(matrix)
    flatmat = get_flatmat(np_matrix, min_duration, pause_rad)
    result = collapse_pause_in_matrix(np_matrix, flatmat)
    return result


def collapse_pause_in_matrix(matrix, flatmat):
    if len(flatmat) == 0:
        return matrix
    else:
        output_mat = []
        if flatmat[0][0] > 1:            
            output_mat = list(matrix[0:flatmat[0][0]-1,])

        for i in range(len(flatmat)):
            start_idx = flatmat[i][0]
            end_idx = flatmat[i][1]            

            result = collapse_to_pause(matrix[start_idx:end_idx+1])
            output_mat.append(result)            

            if i+1 < len(flatmat) and flatmat[i][1] < flatmat[i+1][1] - 1:    
                another_start = flatmat[i][1]+1
                another_end = flatmat[i+1][0]
                rows = matrix[another_start:another_end]
                
                for row in rows:
                    output_mat.append(row.tolist())                    

        if flatmat[-1][1] < len(matrix)-1:    

            the_idx = flatmat[-1][1]+1
            rows = matrix[the_idx:len(matrix),]
            
            for row in rows:                
                output_mat.append(row.tolist())


    if len(output_mat) == 1:
        return output_mat
    else:
        # Group together adjacent pauses, averaging their location
        flatmat2 = []
        output_mat2 = []
        collapse = False
        cstart = 0
        
        for i in range(1, len(output_mat)):
            if output_mat[i][0] != 2 and not collapse:
                continue
            elif output_mat[i][0] != 2 and collapse:
                collapse = False
                flatmat2.append([cstart, i-1])
            elif output_mat[i][0] == 2 and output_mat[i-1][0]==2 and not collapse:                
                cstart = i-1
                collapse = True
            elif output_mat[i][0] == 2 and collapse:
                continue

    if collapse and output_mat[-1][0] == 2:
        flatmat2.append([cstart, len(output_mat)-1])

    if flatmat2 == []:
        output_mat2 = output_mat
    else:        
        output_mat2 = []
        if flatmat2[0][0] > 1:
            end_idx = flatmat2[0][0]            
            output_mat2 = output_mat[0: end_idx]

        for i in range(len(flatmat2)):
            start_idx = flatmat2[i][0]            
            end_idx = flatmat2[i][1]+1
            temp = output_mat[start_idx:end_idx]            
            result = collapse_to_pause(temp)
            output_mat2.append(result)

            if i < len(flatmat2)-1 and flatmat2[i][1] < flatmat2[i+1][0]-1:
                init_idx = flatmat2[i][1]+1
                end_idx = flatmat2[i+1][0]
                result2 = output_mat[init_idx:end_idx]
                temp_result2 = output_mat[init_idx:end_idx+1]
                output_mat2 += result2

        if flatmat2[-1][1] < len(output_mat):    
            begin = (flatmat2[-1][1]+1)
            end = len(output_mat)
            temp = output_mat[begin:len(output_mat)]
            for t in temp:
                output_mat2.append(t)

    # Set flight endpoints equal to pause endpoints
    if output_mat2[0][0]== 1 and output_mat2[1][0] == 2:
        output_mat2[0][4] = output_mat2[1][1]
        output_mat2[0][5] = output_mat2[1][2]

    for i in range(1, len(output_mat2)-1):
        if output_mat2[i][0] == 1 and output_mat2[i-1][0] == 2:            
            output_mat2[i][1] = output_mat2[i-1][1]
            output_mat2[i][2] = output_mat2[i-1][2]
        if output_mat2[i][0] == 1 and output_mat2[i+1][0] == 2:            
            output_mat2[i][4] = output_mat2[i+1][1]
            output_mat2[i][5] = output_mat2[i+1][2]
    
    if output_mat2[-2][0] == 2 and output_mat2[-1][0] == 1:
        output_mat2[-1][1] = output_mat2[-2][1]
        output_mat2[-1][2] = output_mat2[-2][2]

    return output_mat2


def get_flatmat(np_matrix, min_duration=300, pause_rad=60):
    index = 0
    incr = 1    
    t_current = np_matrix[index][3]
    collapse = False
    flatmat = []
    
    counter_here = 0    
    while True:        
        sec_time = np_matrix[index+incr][6]
        if sec_time:
            t_next = sec_time
        else:
            t_next = np_matrix[index+incr][3]

        tmp_again = t_next - t_current
        if str(tmp_again) == "nan":
            tmp_again = min_duration
        if tmp_again >= min_duration:
            tmp = index+incr
            if np_matrix[index+incr][0] == 1:
                #combine
                counter_here += 1    
                left_cols = np_matrix[index:(index+incr+1),1:3]
                right_cols = np_matrix[index:(index+incr+1),4:6]    
                combined = np.concatenate((left_cols, right_cols))
                max_rad = get_max_radius(combined)
            else:                
                left_cols = np_matrix[index:(index+incr+1),1:3]
                max_rad = get_max_radius(left_cols)

            if max_rad > pause_rad and not collapse:
                index += 1                
                sec_time = np_matrix[index][6]
                t_current = np_matrix[index][0]
                incr = 0                
                
            if max_rad > pause_rad and collapse:
                is_else = False
                if np_matrix[index+incr-1, 0] == 4:
                    if index < index+incr - 2:
                        flatmat.append([index, index+incr-2])
                else:
                    end_idx = index + incr - 1
                    if index != end_idx:
                        flatmat.append([index, end_idx])
                    is_else = True

                index += incr

                sec_time = np_matrix[index][6]
                
                t_current = np_matrix[index][0]                
                incr = 0
                collapse = False

            if max_rad <= pause_rad:
                collapse = True

        incr += 1
        if index + incr >= len(np_matrix):
            if max_rad <= pause_rad and t_next - t_current >= min_duration:
                flatmat.append([index, index+incr-1])            
            break

    return flatmat

def collapse_to_pause(matrix):    
    if len(matrix) == 1:
        return matrix[0]

    #check if there is None
    for i, row in enumerate(matrix):
        for j, cell in enumerate(row):
            if not cell:
                matrix[i][j] = float("NaN")
    matrix = np.array(matrix)
    center = np.nanmean(matrix[:,1:3], axis=0)
    
    if str(matrix[-1][6]) != "nan":        
        row = [2, center[0], center[1], matrix[0][3], float("NaN"), float("NaN"), matrix[-1][6]]
    else:        
        row = [2, center[0], center[1], matrix[0][3], float("NaN"), float("NaN"), matrix[-1][3]]

    return row

def get_max_radius(matrix):
    center_of_points = np.nanmean(matrix, axis=0)    
    x_list = matrix[:,0]
    y_list = matrix[:,1]
    dists = np.sqrt((x_list - center_of_points[0])**2 + (y_list - center_of_points[1])**2)
    max_dist = np.nanmax(dists)    
    return max_dist


def initialize_params(matrix):
    np_matrix = np.array(matrix)    
    one = np.where(np_matrix[:,0] == 1)
    two = np.where(np_matrix[:,0] == 2)
    three = np.where(np_matrix[:,0] == 3)
    four = np.where(np_matrix[:,0] == 4)
    
    if len(one) > 1:
        ID1p1 = one[0]+1
        condition1 = one[0][len(one[0])-1]
        if len(one[0]) > 0 and condition1 == len(np_matrix)-1:
            ID1p1 = ID1p1[:len(ID1p1)-1]    
        all_timestamp = np.apply_along_axis(np.mean, 1, np_matrix[:,[3,6]])    
        all_x = np_matrix[:, 1]
        all_y = np_matrix[:, 2]

        #which code 1 is followed by 1
        ind11 = ID1p1[np.where(np_matrix[ID1p1,0] == 1)]    

        #which flight(1) followed by pause (2)
        ind12 = ID1p1[np.where(np_matrix[ID1p1,0] == 2)]

        l1 = len(ind11)
        l2 = len(ind12)

        if (l1 + l2) > 0:
            phatall = l2/(l1+l2)

        if (l1+l2) == 0:
            phatall = len(two)/(len(one) + len(two))


        # ------ flight distances & times & pauses -------
        flight_distances = []
        flight_times = []
        pause_times = []

        for row in np_matrix:
            if row[0] == 1:
                distance = np.sqrt((row[1]-row[4])**2 + (row[2] - row[5])**2)
                flight_distances.append(distance)
                times = row[6] - row[3]
                flight_times.append(times)
            if row[0] == 2:
                pause = row[6] - row[3]
                pause_times.append(pause)

        fxs = np_matrix[one][:,1]
        fys = np_matrix[one][:,2]

        fa = np.zeros(len(one[0]))

        yvals = np_matrix[one][:,5] - np_matrix[one][:,2]
        xvals = np_matrix[one][:,4] - np_matrix[one][:,1]

        IDyg0 = np.where(yvals[:] >=0)[0]
        IDxg0 = np.where(xvals[:] >=0)[0]
        IDyl0 = np.where(yvals[:] < 0)[0]
        IDxl0 = np.where(xvals[:] < 0)[0]

        IDgg = list(set(IDyg0) & set(IDxg0))
        IDlg = list(set(IDyg0) & set(IDxl0))
        IDgl = list(set(IDyl0) & set(IDxg0))
        IDll = list(set(IDyl0) & set(IDxl0))
        IDgg.sort()
        IDlg.sort()
        IDgl.sort()
        IDll.sort()

        fa[IDgg] = np.arctan(list(yvals[IDgg]/xvals[IDgg]))
        fa[IDgl] = np.arctan(list(yvals[IDgl]/xvals[IDgl]))+2*math.pi
        fa[IDlg] = np.arctan(list(yvals[IDlg]/xvals[IDlg])) +math.pi
        fa[IDll] = np.arctan(list(yvals[IDll]/xvals[IDll])) +math.pi

        #flight timestamps
        flight_timestamps = np_matrix[one][:,3]

        pxs = np_matrix[two][:,1] #x
        pys = np_matrix[two][:,2] #y
        pts = np_matrix[two][:,3] #pause timestamps

        the_dict = {'ID1': one, "ID2": two, "ID3": three, "ID4":four, "ID1p1": ID1p1, "allts": all_timestamp, "ind11": ind11, "ind12": ind12, "phatall": phatall, "fd": flight_distances, "ft": flight_times, "fa": fa, "fts": flight_timestamps, "pt": pause_times, "pts": pts, "fxs": fxs, "fys": fys, "pxs": pxs, "pys": pys, "allxs": all_x, "allys": all_y}
        return the_dict
    else:
        print("No mobility flights were found in this participant's location data")
        return None

#simulate_mobility_gaps
#impute the missing gaps hot-tech computation
def simulate_mobility_gaps(matrix, obj, wtype, spread_pars):
    ind11 = obj['ind11']
    ind12 = obj['ind12']
    fd = obj['fd']
    ft = obj['ft']
    fts = obj['fts']
    fa = obj['fa']
    pts = obj['pts']
    allts = obj['allts']
    phatall = obj['phatall']
    fxs = obj['fxs']
    fys = obj['fys']
    pxs = obj['pxs']
    pys = obj['pys']
    allxs = obj['allxs']
    allys = obj['allys']
    pt = obj['pt']
    if len(matrix) == 0:
        return matrix

    f_outmat = []

    for i in tqdm(range(len(matrix)), total=len(matrix)):
        if matrix[i][0] == 1:
            cur_x = matrix[i][4]
            cur_y = matrix[i][5]
            f_outmat.append(matrix[i])

        elif matrix[i][0] <= 3:
            cur_x = matrix[i][1]
            cur_y = matrix[i][2]
            f_outmat.append(matrix[i])

        else: #check here
            if i > 0 and i < len(matrix) - 1:
                varmult = 1
                while True:
                    #----- fw ------                    
                    denominator = varmult*(matrix[i][6]-matrix[i][3])        
                    the_mean = np.mean([matrix[i][3],matrix[i][6]])
                    temp_fw = (fts-the_mean)/denominator
                    temp_pw = (pts-the_mean)/denominator
                    temp_allw = (allts-the_mean)/denominator
                    
                    fw = norm.pdf(list(temp_fw))
                    pw = norm.pdf(list(temp_pw))
                    allw = norm.pdf(list(temp_allw))

                    if len(pts)>0 and len(fts)>0 and np.sum(fw)>0 and np.sum(pw)>0:                        
                        break

                    if len(pts) == 0 and len(fts)>0 and np.sum(fw)>0:
                        break

                    if len(fts) == 0 and len(pts)>0 and np.sum(pw)>0:
                        break
                    
                    varmult = varmult*2    
                s11 = np.nansum(allw[ind11])
                s12 = np.nansum(allw[ind12])
                if (s11 + s12) == 0:
                    phatcur = phatall
                else:
                    phatcur = s12/float(s11+s12)

                if wtype == "LI":
                    new_row = [1, cur_x, cur_y, matrix[i][3], matrix[i+1][1], matrix[i+1][2], matrix[i][6]]
                    f_outmat.append(new_row)
                else:
                    rb_out = random_bridge(x0=cur_x,y0=cur_y,x1=matrix[i+1][1],y1=matrix[i+1][2],t0=matrix[i][3],t1=matrix[i][6],fd=fd,ft=ft,fts=fts,fa=fa,fw=fw,probp=phatcur,pt=pt,pts=pts,pw=pw,allts=allts,allw=allw,ind11=ind11,ind12=ind12,i_ind=i,pxs=pxs,pys=pys,fxs=fxs,fys=fys,allxs=allxs,allys=allys,wtype=wtype,canpause=matrix[i-1][0]==1,spread_pars=spread_pars,niter=100)
                    old_foutmat = f_outmat
                    f_outmat  = f_outmat + rb_out

    return f_outmat

def get_weights(fxs, cur_x, fys, cur_y, varmult, spread_pars, fts, cur_t, t1, t0, weight_type):
    if weight_type == "TL":
        calc = spread_pars[0]*(fts-cur_t)/(varmult*(t1-t0))
        t_score = t.pdf(calc, df=spread_pars[1])
        return t_score
    elif weight_type == "GL":
        pow_xs = (fxs-cur_x)**2
        pow_ys = (fys-cur_y)**2
        sum_x_y = np.array(pow_xs + pow_ys, dtype=np.float64)
        
        distance = np.sqrt(sum_x_y)    
        calc = spread_pars[0]*distance/(50*varmult)
        t_score = t.pdf(calc,df=spread_pars[1])
        return t_score

    elif weight_type == "GLR":
        pow_xs = (fxs-cur_x)**2
        pow_ys = (fys-cur_y)**2
        sum_x_y = np.array(pow_xs + pow_ys, dtype=np.float64)
        
        distance = np.sqrt(sum_x_y)    
        calc = spread_pars[0]*distance/(50*varmult)

        left = t.pdf(calc,df=spread_pars[1])

        temp1 = np.abs((fts-cur_t))%(60*60*24)                    
        temp2 = (60*60*24)-np.abs((fts-cur_t))%(60*60*24)

        min_num = np.minimum(temp1, temp2)
        top = np.array(spread_pars[0]* min_num)
        
        bottom = np.array(varmult*(t1-t0))
        top_div_bottom = list(top/bottom)
        
        right = t.pdf(top_div_bottom,df=spread_pars[1])
        fw = left * right
        return fw

def random_bridge(x0,y0,x1,y1,t0,t1,fd,ft,fts,fa,fw,probp,pt,pts,pw,allts,allw,ind11,ind12,i_ind,pxs,pys,fxs,fys,allxs,allys,wtype,canpause,spread_pars,niter=100):
    success = False
    for i in range(niter):
        outmat = []
        cur_x = x0
        cur_y = y0
        cur_t = t0
        t_arrive = t0
        while True:
            varmult = 1
            while True:
                fw = get_weights(fxs, cur_x, fys, cur_y, varmult, spread_pars, fts, cur_t, t1, t0, weight_type=wtype)
                pw = get_weights(pxs, cur_x, pys, cur_y, varmult, spread_pars, pts, cur_t, t1, t0, weight_type=wtype)
                allw = get_weights(allxs, cur_x, allys, cur_y, varmult, spread_pars, allts, cur_t, t1, t0, weight_type=wtype)
                
                if len(pts)>0 and len(fts)>0 and np.sum(fw)>0 and np.sum(pw)>0:
                    break

                if len(pts) == 0 and len(fts)>0 and np.sum(fw)>0:
                    break

                if len(fts) == 0 and len(pts)>0 and np.sum(pw)>0:
                    break

                varmult = varmult*2 
            
            s11 = np.nansum(allw[ind11])
            s12 = np.nansum(allw[ind12])
            if (s11 + s12) == 0:
                phatcur = probp
            else:
                phatcur = s12/float(s11+s12)
            probp = phatcur
            uniform_sample = np.random.uniform(0, 1, 1)
            if canpause and uniform_sample < probp:
                canpause = False
                normalized_pw = pw/np.sum(pw)
                p_samp = np.random.choice(pt, 1, replace=False, p=normalized_pw)
                if (cur_t + p_samp) < 1:
                    next_line = [2, cur_x, cur_y, cur_t, float("Nan"), float("Nan"), cur_t + p_samp]
                    cur_t = cur_t + p_samp
                    outmat.append(next_line)
                else:
                    break
            else:
                canpause = True
                the_list = list(range(0,len(fa)))
                normalized_fw = fw/np.sum(fw)
                IDsamp = np.random.choice(the_list, 1, replace=False, p=normalized_fw)[0]
                a_samp = fa[IDsamp]
                d_samp = fd[IDsamp]
                t_samp = ft[IDsamp]
                if cur_t + t_samp < t1:
                    next_x = cur_x + np.cos(a_samp)*d_samp
                    next_y = cur_y + np.sin(a_samp) * d_samp

                    if str(next_x) == "nan":
                        break
                    next_line = [1, cur_x, cur_y, cur_t, next_x, next_y, cur_t + t_samp]
                    cur_t = cur_t + t_samp
                    cur_x = next_x
                    cur_y = next_y                    
                    outmat.append(next_line)
                    t_arrive = cur_t
                else:
                    break
        if t_arrive > t0:
            success = True
            break

    if not success:
        return [[1, x0, y0, t0, x1, y1, t1]]
    else:
        cur_x = x0 
        cur_y = y0
        for i in range(len(outmat)):    
            if outmat[i][0] == 1:
                outmat[i][1] = cur_x
                outmat[i][2] = cur_y
                cur_t = outmat[i][6]
                w = (t_arrive-cur_t)/(t_arrive - t0)
                outmat[i][4]=outmat[i][4]*w+x1*(1-w)
                outmat[i][5]=outmat[i][5]*w+y1*(1-w)
                cur_x=outmat[i][4]
                cur_y=outmat[i][5]    


            else:
                cur_t=outmat[i][3]
                outmat[i][1]=cur_x
                outmat[i][2]=cur_y

        outmat[-1][6] = t1
        return outmat

#get significant locations
def sig_locs(mobmat, obj, center_rad=125, timezone="", min_pause_time=600):

    np_matrix = np.array(mobmat, dtype=np.float64)
    obj = initialize_params(np_matrix)
    ID2_from_matrix = obj['ID2'][0]
    if len(ID2_from_matrix) == 0:
        print("No pauses in mobmat within function sig locs")
        return None
    elif len(ID2_from_matrix) == 1:
        x_for_outmat = [np_matrix[ID2_from_matrix[1]][1]]
        y_for_outmat = [np_matrix[ID2_from_matrix[1]][2]]
        timepresent_for_outmat = np.array([0])
        home_for_outmat = np.array([1])
        outmat=[x_for_outmat, y_for_outmat, timepresent_for_outmat, home_for_outmat]        
    else:
        pt_div_min = np.array(obj["pt"])/min_pause_time
        ptred = np.floor(pt_div_min)        
        pause_ids = np.where(ptred > 0)[0]        
        if len(pause_ids) < 2:
            print("No pauses long enough in mobmat within function SigLocs!")
            return None
        
        pmat = []
        for i in range(len(ID2_from_matrix)):
            if ptred[i] > 0:
                row_i = np_matrix[ID2_from_matrix[i]][1:3]
                duplicate_row_i = [row_i] * int(ptred[i])                
                pmat += duplicate_row_i

        kmeansk_v = list(range(2, len(pause_ids)+1))
        previous_fit = []
        for i in range(len(kmeansk_v)):
            kmeansk = kmeansk_v[i]
            kmeans_fit = KMeans(n_clusters=kmeansk, random_state=23).fit(pmat)
            the_distances = pdist(kmeans_fit.cluster_centers_, 'euclidean')
            min_dist = np.min(the_distances)
            if min_dist < center_rad*2 or i == len(kmeansk_v):        
                if i > 0:
                    kmeansk = kmeansk_v[i-1]
                    kmeans_fit = previous_fit
                break            
            previous_fit = kmeans_fit

        cluster_centers = kmeans_fit.cluster_centers_
        num_row_centers = len(cluster_centers)

        outmat = [cluster_centers[:,0], cluster_centers[:,1], np.zeros(num_row_centers), np.zeros(num_row_centers)]

    #Determine time spent at these significant locations
    for i in range(len(ID2_from_matrix)):
        for j in range(num_row_centers):
            euc_dist = np.sqrt((np_matrix[ID2_from_matrix[i]][1] - outmat[0][j])**2 + (np_matrix[ID2_from_matrix[i]][2] - outmat[1][j])**2)
            if euc_dist < center_rad:
                outmat[2][j] = outmat[2][j] + obj['pt'][i]

    #Determine which is home (where is the night spent)    
    for i in range(len(ID2_from_matrix)):
        avg_time = (np_matrix[ID2_from_matrix[i]][6]+np_matrix[ID2_from_matrix[i]][3])/2
        if timezone == "":
            hour_of_day = datetime.datetime.fromtimestamp(avg_time).hour
        else:
            hour_of_day = datetime.datetime.fromtimestamp(avg_time, tz=tz.gettz(timezone)).hour
            
        if hour_of_day >= 21 or hour_of_day < 6:
            for j in range(num_row_centers):
                euc_dist = np.sqrt((np_matrix[ID2_from_matrix[i]][1] - outmat[0][j])**2 + (np_matrix[ID2_from_matrix[i]][2] - outmat[1][j])**2)
                if euc_dist < center_rad:                    
                    outmat[3][j] = outmat[3][j] + obj['pt'][i]
    IDmax = np.argmax(outmat[3])    
    outmat[3][IDmax] = 1
    transposed_outmat = np.array(outmat).transpose()
    zero_ids = np.where(transposed_outmat[:,2] ==0)[0]
    new_outmat = np.delete(transposed_outmat, zero_ids, axis=0)
    sorted_outmat = new_outmat[new_outmat[:,2].argsort()[::-1]]
    
    return sorted_outmat

def home_time(matrix, slout, center_rad):     
    IDhome = np.where(slout[:, 3] == 1)[0]
    
    x_center = slout[IDhome, 0]
    y_center = slout[IDhome, 1]    
    
    tot_time = 0
    for i in range(len(matrix)):
        if matrix[i][0]==2 and np.sqrt((matrix[i][1]-x_center)**2+(matrix[i][2]-y_center)**2)<center_rad: 
            tot_time=tot_time+matrix[i][6]-matrix[i][3]
    result = tot_time/60    
    return result

def distance_traveled(matrix):
    dt = 0
    for i in range(len(matrix)):
        if matrix[i][0]==1:
          temp = np.sqrt((matrix[i][4]-matrix[i][1])**2+(matrix[i][5]-matrix[i][2])**2)
          dt += temp    
    return dt

def radius_of_gyration(matrix, interval): 
    IDskip = np.where(matrix[:,0]==4)[0]
    if len(IDskip) > 0:
        matrix = np.delete(matrix, IDskip, axis=0)

    N = len(matrix)
    w_v = np.zeros(N)
    x_v = np.zeros(N)
    y_v = np.zeros(N)

    for i in range(N):
        if matrix[i][0] == 4:
            continue
        if matrix[i][0] == 3:
            x_v[i] = matrix[i][1]
            y_v[i] = matrix[i][2]
            w_v[i] = interval
        if matrix[i][0] == 1:
            x_v[i] = np.mean([matrix[i][1], matrix[i][4]])
            y_v[i] = np.mean([matrix[i][2], matrix[i][5]])
            w_v[i] = matrix[i][6]-matrix[i][3]
        if matrix[i][0] == 2:
            x_v[i] = matrix[i][1]
            y_v[i] = matrix[i][2]
            w_v[i] = matrix[i][6]-matrix[i][3]

    sum_w_v = np.sum(w_v)
    x_avg = np.sum(w_v*x_v)/sum_w_v
    y_avg = np.sum(w_v*y_v)/sum_w_v
    result = np.sqrt(np.sum(((x_v-x_avg)**2+(y_v-y_avg)**2)*w_v)/sum_w_v)
    
    return result


def max_diameter(matrix):
    IDmv = np.where(matrix[:,0] <=2)[0]
    if len(IDmv) < 2:        
        return 0
    else:
        the_distances = pdist(matrix[IDmv,1:3], 'euclidean')
        max_dist = np.max(the_distances)        
        return max_dist

def max_home_distance(matrix, homex, homey): 
    IDmv = np.where(matrix[:,0] <=2)[0]
    if len(IDmv) == 0:
        return None
    else:
        dfhome = np.zeros(len(IDmv))
        for i in range(len(IDmv)):
            temp = np.sqrt((matrix[IDmv[i]][1]-homex)**2+(matrix[IDmv[i]][2]-homey)**2)
            dfhome[i] = temp
    return np.max(dfhome)

def sig_locs_visited(matrix, slout, center_rad): 
    places_visited = np.zeros(len(slout))
    for i in range(len(matrix)):
        if matrix[i][0] <=3:
            for j in range(len(slout)):
                temp = np.sqrt((slout[j][0]-matrix[i][1])**2 + (slout[j][1]-matrix[i][2])**2)
                if temp < center_rad:
                    places_visited[j] = 1

    result = np.sum(places_visited)
    return result

def avg_flight(matrix, avg_type): 
    #avg_type = "length" or "duration"
    num = 0.0
    total = 0.0
    for i in range(len(matrix)):
        if matrix[i][0] == 1:
            added_number = 0
            if avg_type == "length":
                added_number = np.sqrt((matrix[i][4]-matrix[i][1])**2 + (matrix[i][5]-matrix[i][2])**2)
            else:
                added_number = matrix[i][6] - matrix[i][3]
            total += added_number
            num += 1

    if num == 0:
        return 0
    else:
        return total/num

def std_flight(matrix, std_type): 
    #std_type = "length" or "duration" 
    ID1 = np.where(matrix[:,0]==1)[0]
    if len(ID1) <= 1:
        return 0
    else:
        if std_type == "length":
            temp = ((matrix[ID1,5]-matrix[ID1,2])**2 + (matrix[ID1, 4] - matrix[ID1, 1])**2)
            dist = np.sqrt(list((matrix[ID1,5]-matrix[ID1,2])**2 + (matrix[ID1, 4] - matrix[ID1, 1])**2))
            std_result = np.std(dist)
        else:
            std_result = np.std(matrix[ID1, 6] - matrix[ID1,3])
        return std_result


def prob_pause(matrix): 
    t_pause = 0.0
    t_flight = 0.0
    for i in range(len(matrix)):
        diff = matrix[i][6] - matrix[i][3]
        if matrix[i][0] == 1:
            t_flight += diff
        elif matrix[i][0] == 2:
            t_pause += diff

    result = t_pause/(t_pause+t_flight)
    return result

def sig_loc_entropy(matrix, slout, center_rad): 
    tp = np.zeros(len(slout))
    for i in range(len(matrix)):
        if matrix[i][0] ==2:
            for j in range(len(slout)):
                temp = np.sqrt((slout[j][0]-matrix[i][1])**2 + (slout[j][1]-matrix[i][2])**2)
                if temp < center_rad:
                    tp[j] += matrix[i][6]-matrix[i][3]

    total = 0
    sum_tp = np.sum(tp)
    if sum_tp == 0:
        return 0
    else:
        for i in range(len(slout)):
            p = tp[i]/sum_tp
            if p > 0:
                total = total - p*np.log(p)
    return total

def mins_missing(matrix):
    total = 0.0
    for i in range(len(matrix)):
        if matrix[i][0] == 4:
            total += matrix[i][6] - matrix[i][3]

    return total/60.0

def location_at(matrix, the_time):
    for i in range(len(matrix)):
        if matrix[i][0] <= 2:
            if matrix[i][3] <= the_time and matrix[i][6] >= the_time:
                return matrix[i][1], matrix[i][2]
        elif matrix[i][2] == 3:
            if matrix[i][3] == the_time:
                return matrix[i][1], matrix[i][2]
    return None

def day_dist(i1, i2, mobmat, subset_inds_v,subset_start_time_v,center_rad):
    
    multiplier = np.arange(30*60, 60*60*24, 60*60)
    
    checkpoints = subset_start_time_v[i1] + multiplier    
    checkpoints2 = subset_start_time_v[i2] + multiplier
    matrix1 = mobmat[subset_inds_v[i1]]
    matrix2 = mobmat[subset_inds_v[i2]]
    same_place = np.empty(24)
    same_place[:] = np.nan
    for i in range(24):
        loc1 = location_at(matrix1, checkpoints[i])
        if loc1 == None:
            continue
        abs_diff = np.abs(matrix2[:,3]-checkpoints[i])%(60*60*24)        
        smaller = list(np.where(abs_diff < 30*60)[0])
        greater = list(np.where(abs_diff > (60*60*24)-30*60)[0])
        IDs = smaller + greater
        
        if len(IDs) > 0:
            can_be_zero = False    
            for j in range(len(IDs)):                
                if matrix2[IDs[j]][0] <= 3:                
                    can_be_zero = True
                    temp =  np.sqrt((matrix2[IDs[j]][1]-loc1[0])**2 + (matrix2[IDs[j]][2]-loc1[1])**2)
                    if temp < center_rad:
                        same_place[i] = 1
                
                if can_be_zero and np.isnan(same_place[i]):
                    same_place[i] = 0
        else:
            for j in range(len(matrix2)):                
                if not np.isnan(matrix2[j][3]) and not np.isnan(matrix2[j][6]) and matrix2[j][3] < checkpoints2[i] and matrix2[j][6] > checkpoints2[i]:
                    if matrix2[j][0] !=4 and np.sqrt((matrix2[j][1]-loc1[0])**2 + (matrix2[j][2]-loc1[1])**2) < center_rad:
                        same_place[i] = 1
                    else:
                        same_place[i] = 0
                    break
            
        
    if np.isnan(same_place).all():
        return None
    else:
        return np.nanmean(same_place)


def daily_routine_index(idx, mobmat, subset_inds_v, subset_day_of_week_v, subset_start_time_v, center_rad):    
    submat = mobmat[subset_inds_v[idx]]
    daydist_v = np.ones(len(subset_inds_v))
    daydist_v[:] = np.nan

    for i in range(len(subset_inds_v)):
        if i == idx:
            continue
        
        daydist_v[i] = day_dist(idx, i, mobmat, subset_inds_v, subset_start_time_v, center_rad)

    if len(np.where(daydist_v != -1)[0]) == 0:
        circ_score = None
    else:
        circ_score = np.nanmean(daydist_v)

    if subset_day_of_week_v[idx] == "Saturday" or subset_day_of_week_v[idx] == "Sunday":
        IDcompare = np.concatenate([np.where(subset_day_of_week_v == "Saturday")[0], np.where(subset_day_of_week_v == "Sunday")[0]])
    else:
        IDcompare = np.concatenate([np.where(subset_day_of_week_v == "Monday")[0], np.where(subset_day_of_week_v == "Tuesday")[0], np.where(subset_day_of_week_v == "Wednesday")[0], np.where(subset_day_of_week_v == "Thursday")[0], np.where(subset_day_of_week_v == "Friday")[0]])

    day_d_compare = daydist_v[IDcompare]
    not_nan = day_d_compare[~np.isnan(day_d_compare)]
    if len(IDcompare) == 0 or len(not_nan) == 0:
        week_score = None
    else:        
        week_score = np.nanmean(daydist_v[IDcompare])
        

    return circ_score, week_score


def get_mobility_features(mobmat, obj, mobmatmiss, timezone, center_rad, interval):
    print("Get mobility features...")    
    mobmat = np.array(mobmat)
    mobmatmiss = np.array(mobmatmiss)
    slout = sig_locs(mobmat, obj, center_rad, timezone=timezone) #x, y, timepresent, home
    IDhome = np.where(slout[:,3] ==1)[0]
    if len(IDhome) == 0:
        IDhome = 0
    homex = slout[IDhome,0]
    homey = slout[IDhome,1]

    if timezone != "":
        date = datetime.datetime.fromtimestamp(mobmat[0][3], tz=tz.gettz(timezone))
    else:
        date = datetime.datetime.fromtimestamp(mobmat[0][3])
    current_year = date.year
    current_month = date.month
    current_day = date.day
    current_hour = date.hour
    current_minute = date.minute
    current_second = date.second

    weekday_mapping = {0: "Monday", 1: "Tuesday", 2: "Wednesday", 3: "Thursday", 4: "Friday", 5: "Saturday", 6: "Sunday"}

    subset_inds_v = {}
    subset_day_of_week_v = []
    subset_start_time_v = []
    cur_date = (current_year, current_month, current_day)
    daystr_v = [cur_date]
    day_ind=0
    subset_inds = [0]
    if(len(mobmat)<2):
        outmat = None
        return outmat, slout
    else:
        for i in range(1, len(mobmat)):
            next_date_raw = datetime.datetime.fromtimestamp(mobmat[i][3])
            next_date = (next_date_raw.year, next_date_raw.month, next_date_raw.day)

            if next_date == cur_date and i < len(mobmat)-1:
                if str(mobmat[i][6]) != "nan":
                    end_date_raw = datetime.datetime.fromtimestamp(mobmat[i][6])
                    end_date = (end_date_raw.year, end_date_raw.month, end_date_raw.day)

                subset_inds.append(i)
            else:
                subset_inds_v[day_ind] = subset_inds
                cur_weekday = weekday_mapping[date.weekday()]
                subset_day_of_week_v.append(cur_weekday)
                to_timestamp = datetime.datetime.strptime("{}-{}-{}".format(cur_date[0], cur_date[1], cur_date[2]), "%Y-%m-%d").timestamp()
                subset_start_time_v.append(to_timestamp) 
                
                day_ind += 1
                if mobmat[i-1][0] == 2 and (mobmat[i-1][6]-mobmat[i-1][3])/(60*60*24)>1:
                    ii = 1
                    temp_mid = mobmat[i-1][3]+(60*60*24)*ii
                    mid_date_raw = datetime.datetime.fromtimestamp(temp_mid)
                    mid_date = (mid_date_raw.year, mid_date_raw.month, mid_date_raw.day)
                    while mid_date != next_date:
                        subset_day_of_week_v.append(weekday_mapping[mid_date_raw.weekday()])
                        middate_numeric = datetime.datetime.strptime("{}-{}-{}".format(mid_date[0], mid_date[1], mid_date[2]), "%Y-%m-%d").timestamp()
                        subset_start_time_v.append(middate_numeric)
                        subset_inds_v[day_ind] = [i-1]
                        day_ind += 1
                        if mid_date != daystr_v[-1]:
                            daystr_v.append(mid_date)
                        ii += 1
                        temp_mid = mobmat[i-1][3]+(60*60*24)*ii
                        mid_date_raw = datetime.datetime.fromtimestamp(temp_mid)
                        mid_date = (mid_date_raw.year, mid_date_raw.month, mid_date_raw.day)

                cur_date = next_date
                date = next_date_raw
                subset_inds = [i-1, i] 
                if cur_date != daystr_v[-1] and not (len(daystr_v) == len(subset_inds_v) and len(mobmat) == i+1):
                    daystr_v.append(cur_date)
    
    #### Partition mobmatmiss data into daily subsets: create subsetinds_v and daystr_v.
    cur_date_raw = datetime.datetime.fromtimestamp(mobmatmiss[0][3])
    cur_date_miss = (cur_date_raw.year, cur_date_raw.month, cur_date_raw.day)
    cur_time_miss = (cur_date_raw.hour, cur_date_raw.month, cur_date_raw.second)

    miss_subset_inds_v = {}
    miss_daystr_v = [cur_date_miss]
    miss_dayind = 0
    subset_inds = [0]
    for i in range(1, len(mobmatmiss)):
        next_date_raw = datetime.datetime.fromtimestamp(mobmatmiss[i][3])
        next_date_miss = (next_date_raw.year, next_date_raw.month, next_date_raw.day)

        if next_date_miss == cur_date_miss and i < len(mobmatmiss)-1:
            subset_inds.append(i)
        else:
            cur_date_miss = next_date_miss
            miss_subset_inds_v[miss_dayind] = subset_inds
            if cur_date_miss !=    miss_daystr_v[-1] and not (len(miss_daystr_v) == len(miss_subset_inds_v) and i == len(mobmatmiss)-1):
                miss_daystr_v.append(cur_date_miss)
            subset_inds = [i]
            miss_dayind += 1
            if mobmatmiss[i-1][0] == 2 and (mobmatmiss[i-1][6]-mobmatmiss[i-1][3])/(60*60*24)>1:
                ii = 1
                temp_mid = mobmatmiss[i-1][3]+(60*60*24)*ii
                mid_date_raw = datetime.datetime.fromtimestamp(temp_mid)
                mid_date = (mid_date_raw.year, mid_date_raw.month, mid_date_raw.day)
                while mid_date != next_date_miss:
                    miss_subset_inds_v[miss_dayind] = i-1
                    miss_dayind += 1
                    if mid_date != miss_daystr_v[-1]:
                        miss_daystr_v.append(mid_date)
                    ii += 1
                    temp_mid = mobmatmiss[i-1][3]+(60*60*24)*ii
                    mid_date_raw = datetime.datetime.fromtimestamp(temp_mid)
                    mid_date = (mid_date_raw.year, mid_date_raw.month, mid_date_raw.day)

    ##### intersect mobmat and mobmatmiss to ignore missing data
    IDkeep = []
    IDkeepmiss = []
    for i in range(len(daystr_v)):
        for j in range(len(miss_daystr_v)):
            if daystr_v[i] == miss_daystr_v[j]:
                IDkeep.append(i)
                IDkeepmiss.append(j)
                break

    if len(IDkeep) == 0:
        outmat = None
        return outmat, slout

    daystr_v = np.array(daystr_v)[IDkeep]
    subset_day_of_week_v = np.array(subset_day_of_week_v)[IDkeep]
    subset_start_time_v = np.array(subset_start_time_v)[IDkeep]
    miss_daystr_v = np.array(miss_daystr_v)[IDkeep]

    num_of_features = 15
    outmat = np.zeros((len(daystr_v), num_of_features))
    submat = [] 
    temp_last = mobmat[subset_inds_v[IDkeep[len(daystr_v)-1]]]

    for i in range(len(daystr_v)):
        submat = mobmat[subset_inds_v[IDkeep[i]]]

        if submat[0][0] == 2 and submat[0][3] < subset_start_time_v[i]:
            submat[0][3] = subset_start_time_v[i]

        if submat[-1][0] == 2 and submat[-1][6] > subset_start_time_v[i]+60*60*24:
            submat[-1][6] = subset_start_time_v[i]+60*60*24

        submat_miss = mobmatmiss[miss_subset_inds_v[IDkeep[i]]]
        if submat_miss[0][0] == 4 and submat_miss[0][3] < subset_start_time_v[i]:
            submat_miss[0][3] = subset_start_time_v[i]

        if submat_miss[-1][0] == 4 and submat_miss[-1][6] > subset_start_time_v[i]+60*60*24:
            submat_miss[-1][6] = subset_start_time_v[i]+60*60*24

        if len(submat) == 0 or len(np.where(slout[:,3] ==1)[0]) == 0:
            outmat[i] = [float("NaN")]*15
            outmat[i][12] = 1440
        else:
            outmat[i][0] = home_time(submat, slout, center_rad=200)
            outmat[i][1] = distance_traveled(submat)
            outmat[i][2] = radius_of_gyration(submat, interval)
            outmat[i][3] = max_diameter(submat)
            outmat[i][4] = max_home_distance(submat, homex, homey)
            outmat[i][5] = sig_locs_visited(submat, slout, center_rad)
            outmat[i][6] = avg_flight(submat, "length")
            outmat[i][7] = std_flight(submat, "length")
            outmat[i][8] = avg_flight(submat, "duration")
            outmat[i][9] = std_flight(submat, "duration")
            outmat[i][10] = prob_pause(submat)
            outmat[i][11] = sig_loc_entropy(submat, slout, center_rad)
            outmat[i][12] = mins_missing(submat_miss)
            dri_output = daily_routine_index(i, mobmat, subset_inds_v, subset_day_of_week_v, subset_start_time_v, center_rad)
            outmat[i][13] = dri_output[0] #seven days equally circadian routine
            outmat[i][14] = dri_output[1] #separately


    return outmat, slout, daystr_v #daystr_v = row names
