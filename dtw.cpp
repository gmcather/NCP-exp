#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <pthread.h>
#include <time.h>

using namespace std;

float check(bool condition, float ro, float re){
    return condition ? ro : re;
}

float precision(vector<vector<float> > &sim){
    int cnt = 0;
    for(int i = 0; i < sim.size(); i++){
        int id = -1;    float min_sim = 0.0;
        for(int j = 0; j < sim[i].size(); j++)
            if(min_sim > sim[i][j]){
                min_sim = sim[i][j];
                id = j;
            }
        if(id == i) cnt++;
    }
    return cnt*1.0/sim.size();
}

//m: 18, 35
vector<float> DTW(vector<vector<float> > &seq1, vector<vector<float> > &seq2, int shift, int m=14, int tao=1, float k=0.1, float ro=2, float re=0.5){
    vector<vector<int> > ha(12, vector<int>(12, 0));
    for(int i = 0; i < 12; i++)
        for(int j = 0; j < 12; j++)
            ha[i][j] = (i+j)%12;
    
    int k_size_1 = (int)(seq1.size()*k);
    int k_size_2 = (int)(seq2.size()*k);
    
    vector<vector<float> > cost(seq1.size(), vector<float>(seq2.size(), 0.0));
    vector<vector<float> > prefix_sum(seq1.size(), vector<float>(seq2.size(), 0.0));
    for(int i = 0; i < seq1.size(); i++)
        for(int j = 0; j < seq2.size(); j++){
            if(i>=tao && j>=tao)    prefix_sum[i][j] += prefix_sum[i-tao][j-tao];
            float total_sum = 0;
            for(int k = 0; k < seq1[0].size(); k++)
                total_sum += (seq1[i][k]-seq2[j][ha[k][shift]])*(seq1[i][k]-seq2[j][ha[k][shift]]);
            prefix_sum[i][j] += total_sum;
        }
    
    for(int i = 0; i < seq1.size(); i++)
        for(int j = 0; j < seq2.size(); j++){
            cost[i][j] = prefix_sum[i][j];
            if(i>=m*tao && j>=m*tao)    cost[i][j] -= prefix_sum[i-m*tao][j-m*tao];
            int l1 = i, l2 = j;
            if(l1 > m)  l1 = m;
            if(l2 > m)  l2 = m;
            cost[i][j] /= (l1*l2);
        }
    
    float lmax=0, smax=0, qmax=0;
    vector<float> row(seq1.size());
    vector<float> col(seq2.size());
    for(int i = 0; i < seq1.size(); i++){
        vector<float> tmp(seq2.size(), 0);
        for(int j = 0; j < seq2.size(); j++)
            tmp[j] = cost[i][j];
        nth_element(tmp.begin(), tmp.begin() + k_size_2, tmp.end());
        //sort(tmp.begin(), tmp.end());
        row[i] = tmp[k_size_2];
    }
    for(int j = 0; j < seq2.size(); j++){
        vector<float> tmp(seq1.size(), 0);
        for(int i = 0; i < seq1.size(); i++)
            tmp[i] = cost[i][j];
        nth_element(tmp.begin(), tmp.begin() + k_size_1, tmp.end());
        //sort(tmp.begin(), tmp.end());
        col[j] = tmp[k_size_1];
    }
    
    vector<vector<float> > q(seq1.size(), vector<float>(seq2.size(), 0));
    for(int i=0; i<seq1.size(); i++){
        for(int j=0; j<seq2.size(); j++){
            if(cost[i][j] <= row[i] &&  cost[i][j] <= col[j]){
                for(int di = 1; di <= 2; di++)
                    for(int dj = 1; dj <= 2; dj++){
                        if(di==2 && dj==2)  continue;
                        if(i-di>=0 && j-dj>=0)    q[i][j] = max(q[i-di][j-dj]+1, q[i][j]);
                    }
            }else{
                if(i-1>=0 and j-1>=0)
                    q[i][j] = max(q[i-1][j-1]-check(cost[i-1][j-1]<=row[i-1] && cost[i-1][j-1]<=col[j-1], ro, re), q[i][j]);
                if(i-2>=0 and j-1>=0)
                    q[i][j] = max(q[i-2][j-1]-check(cost[i-2][j-1]<=row[i-2] && cost[i-2][j-1]<=col[j-1], ro, re), q[i][j]);
                if(i-1>=0 and j-2>=0)
                    q[i][j] = max(q[i-1][j-2]-check(cost[i-1][j-2]<=row[i-1] && cost[i-1][j-2]<=col[j-2], ro, re), q[i][j]);
            }
            qmax = max(q[i][j], qmax);
        }
    }
    
    qmax = sqrt(seq2.size())/(qmax>1.0?qmax:1.0);
    
    vector<float> res;
    res.push_back(lmax);
    res.push_back(smax);
    res.push_back(qmax);
    return res;
}

void normalize(vector<vector<float> > &mat){
    for(int i = 0; i < mat.size(); i++){
        float max_col_num = 0.0;
        for(int j = 0; j < mat[i].size(); j++)
            max_col_num = max(max_col_num, mat[i][j]);
        for(int j = 0; j < mat[i].size(); j++)
            if(max_col_num>0)   mat[i][j] /= max_col_num;
    }
}

vector<vector<float> > read_from_csv(string infile){
    char str[200000];
    vector<vector<float> > mat;
    int row_size = 0;
    ifstream ifs(infile.c_str());
    while(ifs.getline(str, 100000)){
        string sstr(str);
        stringstream ss(sstr);
        vector<float> row;
        char num[50];
        while(ss.getline(num, 50, ','))
            row.push_back(atof(num));
        
        if(row_size==0) mat.resize(row.size());
        for(int i = 0; i < row.size(); i++)
            mat[i].push_back(row[i]);
        row_size++;
    }
    normalize(mat);
    return mat;
}

vector<int> opt_shift(vector<vector<float> > &seq1, vector<vector<float> > &seq2){
    vector<float> global_1(12, 0), global_2(12, 0);
    for(int i = 0; i < seq1.size(); i++)
        for(int j = 0; j < 12; j++)
            global_1[j] += seq1[i][j];
    for(int i = 0; i < seq2.size(); i++)
        for(int j = 0; j < 12; j++)
            global_2[j] += seq2[i][j];
    
    vector<pair<float, int> > shifts;
    for(int i = 0; i < 12; i++){
        float tmp = 0;
        for(int k = 0; k < 12; k++)
            tmp += global_1[k]*global_2[(i+k)%12];
        shifts.push_back(make_pair<float, int>(tmp, i));
    }
    
    for(int i = 0; i < 12; i++)
        for(int j = i+1; j < 12; j++)
            if(shifts[i].first < shifts[j].first)
                swap(shifts[i], shifts[j]);
    
    vector<int> res;
    for(int i = 0; i < 3; i++)
      res.push_back(shifts[i].second);
    return res;
}

int main(int argc, const char * argv[]) {
    vector<vector<float> > seq1 = read_from_csv(argv[1]);
    vector<vector<float> > seq2 = read_from_csv(argv[2]);
    vector<int> shifts = opt_shift(seq1, seq2);
    vector<float> ret = DTW(seq1, seq2, shifts[0], atoi(argv[3]));
    cout << ret[2] << endl;
    return 0;
}
