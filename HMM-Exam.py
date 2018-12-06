
# coding: utf-8

# In[18]:


from Bio import SeqIO
import random
import math
import copy


# In[19]:


class HMM():
    def __init__(self,one_one, one_two, two_one, two_two):
        self.num_states = 2
        self.prior = [0.996 , 0.004]
        self.transition = [[one_one, one_two], [two_one, two_two]]
        self.emission = [{"A": 0.291, "T": 0.291, "C": 0.209, "G": 0.209},
                         {"A": 0.169, "T": 0.169, "C": 0.331, "G": 0.331}]
        self.path=[]

    def logprob(self, sequence, states):
        probability =[]
        current_state = states[0]
        current_obs = sequence[0]
        probability = math.log(self.prior[current_state]) + math.log(self.emission[current_state][current_obs])
        for index in range(1,len(sequence)):
            current_obs = sequence[index]
            current_state = states[index]
            prev_state = states[index -1]
            trans_prob = math.log(self.transition[prev_state][current_state])
            emission_prob = math.log(self.emission[current_state][current_obs])
            probability =  trans_prob + emission_prob + probability
        return probability

    def viterbi(self, sequence):
        steps = [[0 for x in range(2)] for y in range(len(sequence))]
        prev_table=  [[0 for x in range(2)] for y in range(len(sequence))]
        prev_table[0][0] = -1
        prev_table[0][1] = -1

        steps[0][0]= math.log(self.prior[0]) + math.log(self.emission[0][sequence[0]])
        steps[0][1] = math.log(self.prior[1]) + math.log(self.emission[1][sequence[0]])


        for sym in range(1,len(sequence)):
            emission_H = math.log(self.emission[1][sequence[sym]])
            emission_L = math.log(self.emission[0][sequence[sym]])
            L_to_L = self.prob_of_path((0,0),
                                       steps[sym-1][0],
                                       sequence[sym]
                                       )
            H_to_L = self.prob_of_path((1,0),
                                       steps[sym-1][1],
                                       sequence[sym])
            H_to_H = self.prob_of_path((1,1),
                                       steps[sym-1][1],
                                       sequence[sym])
            L_to_H = self.prob_of_path((0,1),
                                       steps[sym-1][0],
                                       sequence[sym])

            steps[sym][0] = emission_L + max(L_to_L , H_to_L)
            steps[sym][1] = emission_H + max(L_to_H , H_to_H)
            prev_table[sym][0] =  0 if  max(L_to_L, H_to_L) == L_to_L else 1
            prev_table[sym][1] =  0 if max(L_to_H , H_to_H) == L_to_H else 1
        return self.back_track(steps,prev_table)

    def back_track(self,prob_table,prev_table):
        path = []
        last_sym = prob_table[len(prob_table)-1]
        start_index = 0
        if(max(last_sym[0],last_sym[1])) == last_sym[0]:
            start_index = 0
            path.append(0)
        else:
            start_index = 1
            path.append(1)
        for index in range(len(prev_table)-1,0,-1):

            prev_node = prev_table[index][start_index]
            path.append( prev_node)
            start_index = prev_node
        path.reverse()
        return path

    def prob_of_path(self,path,prev_probablity,sym):
        start = path[0]
        end = path[1]
        return prev_probablity + math.log(self.transition[start][end])

def transition_prob(sequence):
    prev = 0
    one_one, one_two, two_one, two_two = 0,0,0,0
    for i in range(len(sequence)):
        if sequence[i]==0 and prev==0:
            one_one +=1
        elif sequence[i]==1 and prev==0:
            two_one +=1
        elif sequence[i]==0 and prev==1:
            one_two +=1
        elif sequence[i]==1 and prev==1:
            two_two +=1
        prev = sequence[i]      
    return one_one, one_two, two_one, two_two

def states_and_segment(sequence):
    AT, GC, AT_seg, GC_seg = 0,0,0,0
    prev = -1
    if sequence[1]==0:
        AT_seg +=1
    else:
        GC_seg +=1
    for i in range(len(sequence)):
        if sequence[i]==0:
            AT +=1
        else:
            GC +=1
        if sequence[i]==1 and prev==0:
            GC_seg +=1
        elif sequence[i]==0 and prev==1:
            AT_seg +=1
        prev = sequence[i]
    return AT, GC, AT_seg, GC_seg


# In[20]:


#Driver Code
string = ""
fasta_sequences = SeqIO.parse(open(r"C:\Users\Chandrima\Downloads\sequence.fasta"),'fasta')
for fasta in fasta_sequences:
        name, sequence = fasta, str(fasta.seq)
string = str(name) + "\n"
curr = "\n" +"Transition probabilities" + "\n" + "One_One Two_One One_Two Two_Two" + "\n"
oo,ot, to,tt = 0.999,0.001,0.01,0.99
string += curr
length = viterbi.__len__()
i=10
states, segments = "", ""
while i!=0:
    hmm = HMM(oo,ot, to,tt)
    viterbi = hmm.viterbi(sequence)
    logprob = hmm.logprob(sequence, viterbi)
    one_one, one_two, two_one, two_two = transition_prob(viterbi)
    AT, GC, AT_seg, GC_seg = states_and_segment(viterbi)
    oo,ot, to,tt = one_one/(one_one+one_two), one_two/(one_one+one_two), two_one/(two_one+two_two), two_two/(two_one+two_two)
    oo,ot, to,tt = round(oo,6) ,round(ot,6) , round(to,6), round(tt,6)
    curr = str(oo) + " " + str(ot) + " " + str(to) + " "+ str(tt) + "\n"
    states += str(AT) + " " + str(GC) + "\n"
    segments += str(AT_seg) + " " + str(GC_seg) + "\n"
    string += curr
    i -=1

string = string + "\n" + "Final States Calculated after each iteration" + "\n" + "AT rich" + " " + "GC rich" + "\n" + states 
string = string + "\n" + "Final Segments Calculated after each iteration" + "\n" + "AT rich" + " " + "GC rich" + "\n" + segments 
string = string + "\n" + "GC Rich Region" + "\n"
prev = 0
start_loc, end_loc = -1, -1
for i in range(len(viterbi)):
    if viterbi[i]==1 and prev==0:
        start_loc = i+1
    if viterbi[i]==0 and (prev==1 or i==len(viterbi)-1):
        end_loc = i
    if start_loc>=0 and end_loc>0:
        curr =  str(start_loc) + " "+ str(end_loc) + "\n"
        string = string + curr
        start_loc, end_loc = -1, -1
    prev = viterbi[i]
print (string)
text_file = open("HMM_output_exam.txt", "w")
text_file.write(string)
text_file.close()

