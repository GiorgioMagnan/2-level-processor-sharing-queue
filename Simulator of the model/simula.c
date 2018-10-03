#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<assert.h>
#include "RngStream.h"

#define MAXQUEUESIZE 5000

/*Arrival rate*/
#define LAMBDA 0.8

/*Service rate*/
#define MU 1


/*Job size*/
#define MEANSIZE 1

#define EXPONENTIAL 1
#define UNIFORM 2
#define BIMODAL 3

typedef enum {ARRIVAL, DEPARTURE} event_type; 


RngStream g;

double simtime;
double queue1[MAXQUEUESIZE];
double queue2[MAXQUEUESIZE];
double arrivaltimes[MAXQUEUESIZE];
double alpha2vec[MAXQUEUESIZE];
double endqueue1[MAXQUEUESIZE];

int first[MAXQUEUESIZE];
int tagged[MAXQUEUESIZE];

unsigned int jobs;
unsigned int jobs1;
unsigned int jobs2;
double next_arrival;

int size_dist=BIMODAL;

double th = (MEANSIZE*6.0);

/*Statistics*/
double totjob;
unsigned completions;
double totresp;

double nbatches;
double totbatches;
double maxbatchsize;

double totbatchtagged;
double nbatchtagged;

double currentbusy;
double totalbusy;

double empty1;                  /*compute the total idle time of Q1*/
unsigned long totalarrivals;    /*total number of arrived jobs*/
unsigned long biggerarrivals;   /*computed the expected number of jobs entering in Q2*/
double arrivalsize;             /*compute the expected size of the arrivals*/
double startempty;              /*beginning of latest idle period in Q1*/

double alpha_one;
double alpha_two;
double trsmaller;
double trbiggerq2;

double samplewaiting;
unsigned long sampleswn;


double exponential(double rate) {
    return -log(RngStream_RandU01(g))/rate;
}

double uniform(double min, double max) {
    return RngStream_RandU01(g)*(max-min)+min;
}

double bimodal(double min, double max, double p) {
    return (RngStream_RandU01(g)<p) ? exponential(1/min) : exponential(1/max);
}

double samplesize(double mean) {
    switch(size_dist) {
        case EXPONENTIAL:
            return exponential(1/mean);
        case UNIFORM:
            return uniform(0,1*mean);
        case BIMODAL:
            return bimodal(0.2, 1.8, 0.5);
    }
    return mean;
}


void initialize() {
    unsigned long seed[6];
    int i;
    g = RngStream_CreateStream("g1");
    simtime = 0.0;
    
    do {
        printf("Setting seed...");
        seed[0] = time(NULL) % 4294944443;
        seed[1] = seed[0] +1;
        seed[2] = seed[1] +1;
        seed[3]=seed[4]=seed[5]=seed[0];
    } while(RngStream_SetSeed (g, seed)<0);
    
    for(i =0; i < MAXQUEUESIZE; i++) {
        first[i] = 0;
        tagged[i] = 0;
        alpha2vec[i] = 0.0;
        endqueue1[i] = 0.0;
    }

    jobs = 0;
    jobs1 = 0;
    jobs2 = 0;
    next_arrival = exponential(LAMBDA);
    
    totjob = 0.0;
    totresp = 0.0;
    completions = 0;
    
    nbatches = 0.0;
    totbatches = 0.0;
    totbatchtagged = 0.0;
    nbatchtagged = 0.0;
    maxbatchsize = 0.0;
    
    currentbusy = 0.0;
    totalbusy = 0.0;
    
    empty1 = 0.0;
    
    arrivalsize = 0.0;
    totalarrivals = 0;
    biggerarrivals = 0;
    startempty = 0.0;
    
    alpha_one = 0.0;
    alpha_two = 0.0;
    
    samplewaiting = 0.0;
    sampleswn = 0;
    trsmaller = 0.0;
    trbiggerq2 = 0.0;
    
    printf("  done\n");
}


double next_event(event_type *type) {
    double next_event = next_arrival;
    *type = ARRIVAL;
    double completion;
    
    int i, position=0;
    
    if (jobs1 > 0) {
        for (i=0; i<jobs; i++) {
/*            assert(queue1[i]>=0.0);
            assert(queue2[i]>=0.0); */
            completion = queue1[i] * jobs1 / MU;
            if (queue1[i] > 0) {
                if (completion + simtime < next_event) {
                    next_event = completion + simtime;
                    *type = DEPARTURE;
                    position = i;
                }
            }
        }
    }
    else if (jobs2 > 0) {
        for (i=0; i<jobs; i++) {
            completion = queue2[i] * jobs2 / MU;
            if (completion + simtime < next_event && queue1[i]==0) {
                next_event = completion + simtime;
                *type = DEPARTURE;
                position = i;
            }
        }
    }
    
    
    
    if (*type == DEPARTURE) {
        double flip;
        int flip2;
        
        flip = queue1[position];
        queue1[position]=queue1[jobs-1];
        queue1[jobs-1]=flip;
        
        flip = queue2[position];
        queue2[position]=queue2[jobs-1];
        queue2[jobs-1]=flip;
        
        flip = arrivaltimes[position];
        arrivaltimes[position] = arrivaltimes[jobs-1];
        arrivaltimes[jobs-1] = flip;
        
        flip = endqueue1[position];
        endqueue1[position] = endqueue1[jobs-1];
        endqueue1[jobs-1] = flip;
        
        flip2 = first[position];
        first[position] = first[jobs-1];
        first[jobs-1] = flip2;

        flip2 = tagged[position];
        tagged[position] = tagged[jobs-1];
        tagged[jobs-1] = flip2;

    }
    
    return next_event;
}


void process_event(event_type t, double ts) {
    int i;
    double elapsed = ts - simtime;
    double workdone = (jobs1 > 0) ? elapsed * MU / jobs1 : elapsed * MU /jobs2;
    double size;
    
    totjob += jobs*elapsed;
    
    switch (t) {
        case ARRIVAL:
            assert(jobs<MAXQUEUESIZE);
            assert(jobs==jobs1+jobs2);
            
            totalarrivals++;
            
            if (jobs1 == 0) {
                currentbusy = ts;
                empty1 += ts-startempty;
            }
            
            /*update residual job sizes*/
            if (jobs1 > 0) /* Queue1 is working */
                for (i=0; i<jobs; i++)
                    queue1[i] -= (queue1[i] > 0.0 )? workdone : 0.0;
            else
                if (jobs2 > 0) { /* Queue2 is working */
                    /*update residual work*/
                    for (i=0; i<jobs; i++) {
                        assert(queue1[i]==0.0);
                        queue2[i] -= (queue2[i] > 0.0)? workdone : 0.0;
                    }
                    
                    /*update alpha2*/
                    for (i=0; i<jobs; i++) {
                        alpha2vec[i] += (queue2[i] > 0 && queue1[i] == 0) ? ts-simtime : 0;
                    }
                }
            
            if (totalarrivals%1000 == 0) {/*sample the waiting time for the job*/
                double tottimeq1 = 0.0;
                for (i=0; i<jobs; i++) {
                    tottimeq1+=queue1[i];
                }
                
                samplewaiting += tottimeq1;
                sampleswn++;
            }
            
            /*Add the job in queue*/
            size = samplesize(MEANSIZE);
            
            alpha2vec[jobs] = 0.0;


            if (size > th) {
                queue1[jobs] = th;
                queue2[jobs] = size - th;
                arrivalsize += th;
                biggerarrivals++;
                first[jobs] = 1;
                tagged[jobs] = uniform(0, 1) > 0.99 ? 1 : 0;
            }
            else {
                queue1[jobs] = size;
                queue2[jobs] = 0.0;
                arrivalsize += size;
                first[jobs] = 0;
            }
            arrivaltimes[jobs++] = ts;
            if (th > 0.0)
                jobs1++;
            else
                jobs2++;
            next_arrival = ts + exponential(LAMBDA);
            
            break;
        case DEPARTURE:
            if (jobs1 == 1) { /*end of busy period of queue1*/
                double currentbatch = 0.0;
                int foundtagged = 0;

                nbatches += 1;
                totalbusy += ts - currentbusy;
                startempty = ts;


                //ALPHAONE
                for(i = 0; i < jobs; i++) {
                    if(first[i] == 1) {
                        alpha_one += ts - arrivaltimes[i];
                        currentbatch  += 1.0;
                        first[i] = 0;
                        endqueue1[i] = ts;
                        if (tagged[i]) {
                            foundtagged = 1;
                            tagged[i] = 0;
                        }
                    }
                }
                    
                if (currentbatch > maxbatchsize)
                    maxbatchsize = currentbatch;
                
                totbatches += currentbatch;
                
                if (foundtagged) {
                    totbatchtagged += currentbatch;
                    nbatchtagged += 1;
                }
            }
            else if (jobs1 == 0 && jobs2 > 0){
                /*update alpha_two*/
                for (i=0; i<jobs; i++)
                    alpha2vec[i] += (queue2[i] > 0 && queue1[i] == 0) ? ts-simtime : 0;
            }
            
            /*Completion at queue1 and departure */
            if (queue1[jobs-1]>0.0 && queue2[jobs-1]==0.0 && jobs1>0) {
                jobs--;
                jobs1--;
                completions++;
                totresp+= (ts - arrivaltimes[jobs]);
                trsmaller += (ts - arrivaltimes[jobs]);
                
                /*update residual job sizes in queue1*/
                for (i=0; i<jobs; i++)
                    queue1[i] -= (queue1[i]>0) ? workdone : 0.0;
                
            }
            /*Completion at queue1 and move to queue2*/
            else if (queue1[jobs-1]>0 && queue2[jobs-1]>0 && jobs1>0) {
                queue1[jobs-1] = 0.0;
                jobs1--;
                jobs2++;
    
                /*update residual job sizes in queue1*/
                for (i=0; i<jobs; i++)
                    queue1[i] -= (queue1[i]>0) ? workdone : 0.0;
             


            }
            /*Completion at queue2*/
            else if (queue1[jobs-1]==0 && queue2[jobs-1]>0 && jobs1==0 && jobs2>0) {
                jobs--;
                jobs2--;
                completions++;
                totresp += (ts - arrivaltimes[jobs]);
                
                /*update residual job sizes in queue2*/
                for (i=0; i<jobs; i++)
                    queue2[i] -= (queue2[i]>0 && queue1[i]==0.0) ? workdone : 0.0;
                queue2[jobs] = 0.0;
                
                /*update alpha_two*/
                alpha_two += alpha2vec[jobs];
                alpha2vec[jobs] = 0.0;
                
                /*Rescaled alpha2_two*/
                trbiggerq2 += ts - endqueue1[jobs];
                endqueue1[jobs] = 0.0;
            }
            /*This is an error*/
            else {
                printf("\nError: something is going wrong\n");
                exit(1);
            }
        
            
    }
    
    simtime = ts;
}


void print_state() {
    int i;
    printf("\n--------------------\n");
    printf("Simtime: %f \n", simtime);
    printf("Next arrival: %f\n", next_arrival);
    for (i=0; i<jobs; i++)
        printf("%d: %f;", i+1, queue1[i]);
    
    printf("\n--------------------\n");
}


int main() {
    initialize();
    unsigned long events = 0;
    
    double ne;
    event_type e;
    
    for (events=0; events < 90000000; events++) {
        /*	print_state();*/
        ne = next_event(&e);
        process_event(e, ne);
    }
    
    printf("\nAverage customers (N): %f", totjob/simtime);
    printf("\nExpected response time (R): %f\n", totresp/completions);
    printf("\nExpected Batch size at Q2 (a): %f", totbatches/nbatches);
    printf("\nExpected busy period at Q1: %f", totalbusy/nbatches);
    printf("\nExpected empty period at Q1: %f", empty1/nbatches);
    printf("\nExpected size of the jobs in Q1: %f", arrivalsize/totalarrivals);
    printf("\nProbability of an arrival at Q2: %f", biggerarrivals*1.0/totalarrivals);
    printf("\nExpected waiting time at Q1 under FIFO (W1): %f", samplewaiting/sampleswn);
    printf("\nExpected batch size for tagged customer (b+1): %f", totbatchtagged / nbatchtagged);
    printf("\nMaximum batch size: %f", maxbatchsize);
    printf("\nLoad factor at Q1 (Rhoa1): %f\n", 1.0 - empty1/simtime);
    printf("\nAlpha one at Q1: %f", alpha_one / biggerarrivals);
    printf("\nAlpha two at Q2: %f", alpha_two / biggerarrivals);
    printf("\nExpected response time at the second queue for larger jobs: %f", trbiggerq2/biggerarrivals);
    printf("\nExpected response time at the second queue for larger jobs (cross): %f", alpha_two/biggerarrivals * simtime /empty1);
    printf("\nExpected response time for smaller jobs: %f", trsmaller/(totalarrivals-biggerarrivals));
    printf("\nProbability of bigger arrivals: %f\n\n", biggerarrivals*1.0/totalarrivals);
    return 0;
}









