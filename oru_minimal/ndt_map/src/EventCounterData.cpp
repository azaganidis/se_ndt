/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, AASS Research Center, Orebro University.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef CEVENT_COUNTER_DATA_H
#define CEVENT_COUNTER_DATA_H
#include<cmath>
#include <Eigen/Dense>
#include <inttypes.h>

///Defines for presets
#define EVENTMAP_OCCU  255
#define EVENTMAP_UNKNOWN  127
#define EVENTMAP_FREE  0
#define EVENTMAP_NUMBER_OF_EVENTS_STORED 128


///Recency filtering stuff
#define EVENTMAP_USE_RECENCY_FILTERING
#define EVENTMAP_OBSERVATION_LIMIT 30000.0
#define Wforget (5000.0/5001.0)

/**
* Structure in the map
**/
struct TEventData
{
    uint8_t 	occval; 			///<Occupancy value
    float 		a_exit_event; ///Number of exit events ( OCC 2 EMP)
    float  		b_exit_event; ///Number of times the cell is perceived as occupied

    float 		a_entry_event;
    float  		b_entry_event; ///Number of times the cell is perceived as empty
    uint64_t  events;        /// A storage for the last 64 events (as bits)

    TEventData():
        occval(127),
        a_exit_event(1.0),
        b_exit_event(1.0),
        a_entry_event(1.0),
        b_entry_event(1.0),
        events(0)
    {
    };

    ~TEventData()
    {
    }


    const TEventData &operator =(const TEventData &copy)
    {
        occval = copy.occval;
        a_exit_event = copy.a_exit_event;
        b_exit_event = copy.b_exit_event;
        a_entry_event = copy.a_entry_event;
        b_entry_event = copy.b_entry_event;
        events = copy.events;

        return *this;
    }

    TEventData(const TEventData& copy):
        occval(copy.occval),
        a_exit_event(copy.a_exit_event),
        b_exit_event(copy.b_exit_event),
        a_entry_event(copy.a_entry_event),
        b_entry_event(copy.b_entry_event),
        events( copy.events)
    {
    }


    void Reset()
    {
        occval = 127;
        a_exit_event = 1;
        b_exit_event = 1;
        a_entry_event = 1;
        b_entry_event = 1;
    }

    bool getBit(int ind)
    {
        if((events & (1<<ind))>0) return true;
        return false;
    }


    float computeShortTermOccupancy()
    {
        static float OCC_LOG_PROP = log(0.6/0.4);
        static float EMP_LOG_PROP = log(0.3/0.7);
        uint bs=0;

        if( (b_entry_event+b_exit_event-2) > 63) bs=64;
        else bs = b_entry_event+b_exit_event-2;
        uint cnt=0;
        float log_ogg_prob  = 0 ;
        for(uint i = 0; i<bs; i++)
        {
            if(getBit(i))
            {
                log_ogg_prob += OCC_LOG_PROP;
                cnt++;
            }
            else
            {
                log_ogg_prob += EMP_LOG_PROP;
            }
        }
        return (1.0- 1.0/(1.0+exp(log_ogg_prob)));
    }

    int getObservations()
    {
        return (b_entry_event+b_exit_event);
    }

    /**
    * Updates a new measurement
    */
    void updateSimple(uint8_t meas)
    {
        events = events<<1;
        if(meas > EVENTMAP_UNKNOWN)  ///OCC
        {
            events |= 1;
        }
        if(occval > EVENTMAP_UNKNOWN)  ///cell is considered Occupied
        {
            updateExitEvent(meas, 1.0);
        }
        else if(occval < EVENTMAP_UNKNOWN)
        {
            updateEntryEvent(meas,1.0);
        }
        else
        {
            occval = meas; ///Just update
        }
    }

    /**
    * Updates a new measurement
    */
    void updateSimple(uint8_t meas, float w)
    {
        events = events<<1;
        if(meas > EVENTMAP_UNKNOWN)  ///OCC
        {
            events |= 1;
        }
        if(occval > EVENTMAP_UNKNOWN)  ///cell is considered Occupied
        {
            updateExitEvent(meas, w);
        }
        else if(occval < EVENTMAP_UNKNOWN)
        {
            updateEntryEvent(meas,w);
        }
        else
        {
            occval = meas; ///Just update
        }
#ifdef EVENTMAP_USE_RECENCY_FILTERING
        performRecencyFiltering();
#endif
    }

    void performRecencyFiltering()
    {
        if(occval > EVENTMAP_UNKNOWN)  ///We are in Exit state
        {
            if(b_exit_event>EVENTMAP_OBSERVATION_LIMIT)
            {
                float w = EVENTMAP_OBSERVATION_LIMIT / b_exit_event;
                b_exit_event *= w;
                a_exit_event *= w;
            }
            a_entry_event = 1.0 + (a_entry_event-1.0) * Wforget;
            b_entry_event = 1.0 + (b_entry_event-1.0) * Wforget;
        }
        else
        {
            if(b_entry_event>EVENTMAP_OBSERVATION_LIMIT)
            {
                float w = EVENTMAP_OBSERVATION_LIMIT / b_entry_event;
                b_entry_event *= w;
                a_entry_event *= w;
            }
            a_exit_event = 1.0 + (a_exit_event-1.0) * Wforget;
            b_exit_event = 1.0 + (b_exit_event-1.0) * Wforget;
        }
    }


    /**
    * Update Exit event (cell is OCC)
    */
    void updateExitEvent(uint8_t meas, float w)
    {
        if(meas>EVENTMAP_UNKNOWN)  ///OCC
        {
            b_exit_event+=1.0;
        }
        else
        {
            if(w>1.0) fprintf(stderr,"EXIT W larger that one = %f\n ",w);
            a_exit_event+=w;
            b_exit_event+=1.0;
            occval = EVENTMAP_FREE;
        }
    }
    /**
    * Update Entry event (Cell is Free)
    */
    void updateEntryEvent(uint8_t meas, float w)
    {
        if(meas>EVENTMAP_UNKNOWN)  ///OCC
        {
            b_entry_event+=1.0;
            if(w>1.0) fprintf(stderr,"W larger that one = %f\n ",w);
            a_entry_event+=w;
            occval = EVENTMAP_OCCU;
        }
        else
        {
            b_entry_event+=1.0;
        }
    }

    ///Expected number of observations/event
    int getEntryN()
    {
        return (int)((float)b_entry_event/(float)a_entry_event +0.5);
    }
    ///Expected number of observations/event
    int getExitN()
    {
        return (int)((float)b_exit_event/(float)a_exit_event +0.5);
    }


    double entryL()
    {
        if(b_entry_event == 0) fprintf(stderr,"B Zero, which is not possible!!\n");
        return ((double)a_entry_event/(double)b_entry_event );
    }
    double exitL()
    {
        if(b_exit_event == 0) fprintf(stderr,"B Zero, which is not possible!!\n");
        return ((double)a_exit_event/(double)b_exit_event);
    }
    double L()
    {
        if(b_exit_event == 0) fprintf(stderr,"B Zero, which is not possible!!\n");
        return ((double)(a_exit_event+a_entry_event)/(double)(b_exit_event+b_entry_event));
    }

    double fac(int n)
    {
        double t=1;
        for (int i=n; i>1; i--)
            t*=i;
        return t;
    }

    double Bin(int n,double p,int r)
    {
        return fac(n)/(fac(n-r)*fac(r))*pow(p,r)*pow(1-p,n-r);
    }

    float binaryBayesUpdate(float Pold, float P)
    {
        return (( Pold*P) / ( Pold*P + (1.0f-Pold)*(1.0f-P))); 	///< Updated probability
    }

    /**
    * To prevent numerical difficulties
    */
    void normalizeProb(float &p)
    {
        if(p>0.999) p = 0.999;
        if(p<0.001) p = 0.001;

    }
    /////////////////////////////////////////////////////////////////////////////////////////////77
    /// Compute different probabilities for the cell
    /////////////////////////////////////////////////////////////////////////////////////////////77
    /////////////////////////////////////////////////////////////////////////////////////////////77
    /**
    * Compute the probability of this cell being Static occupied.
    */
    float getOccStaticLikelihood()
    {
        if( (b_entry_event + b_exit_event)<20.0 )
        {
            return 0.5f;
        }
        Eigen::Matrix2f P;
        Eigen::Vector2f u1(0.5, 0.5);
        Eigen::Vector2f P0;
        float Lex = exitL();
        float Len = entryL();
        P(0,0) = (1.0-Len);
        P(0,1) =  Len;
        P(1,0)	 = Lex;
        P(1,1) = (1-Lex);

        P0 = u1.transpose() * P;
        return (P0(1));
    }

    /**
    * Returns the probability of static empty
    */
    float getFreeStaticLikelihood()
    {
        if( (b_entry_event + b_exit_event)<20.0 )
        {
            return 0.5f;
        }
        Eigen::Matrix2f P;
        Eigen::Vector2f P0;
        Eigen::Vector2f u1(0.5, 0.5);
        float Lex = exitL();
        float Len = entryL();
        P(0,0) = (1.0-Len);
        P(0,1) =  Len;
        P(1,0)	 = Lex;
        P(1,1) = (1-Lex);

        P0 = u1.transpose() * P;
        return (P0(0));
    }

    /**
    * Returns the probability of the cell having dynamics after 2^N
    * Observations. The larger the N, and higher the value, the more "Semi-static"
    * is the cell. Note that the filter returns
    * 0.5 >= P <=1.0 and P = 0.0 only if there is not enough evidence about the
    * Behaviour
    */
    float computeSemiStaticLikelihood(int N)
    {
        if( (b_entry_event + b_exit_event)<20.0 )
        {
            return 0.0f;
        }


        Eigen::Matrix2f P;
        Eigen::Vector2f u2(0, 1.0);
        Eigen::Vector2f u3(1.0, 0);
        float Lex = exitL();
        float Len = entryL();
        Eigen::Vector2f P2,P3;
        P(0,0) = (1.0-Len);
        P(0,1) =  Len;
        P(1,0)	 = Lex;
        P(1,1) = (1-Lex);
        for(int i=0; i<N; i++) P = P*P;
        P2 = u2.transpose() * P;
        P3 = u3.transpose() * P;

        float Po = P2(1);
        float Pu = binaryBayesUpdate(Po, P3(0));
        normalizeProb(Pu);

        return Pu;
    }

    float getOccupancyNow()
    {
        if( (b_entry_event + b_exit_event)<50.0)
        {
            return 0.5f;
        }
        float Po = 0;


        Po = computeShortTermOccupancy();
        float Lex = exitL();
        float Len = entryL();
        Eigen::Matrix2f P;
        P(0,0) = (1.0-Len);
        P(0,1) =  Len;
        P(1,0)	 = Lex;
        P(1,1) = (1-Lex);
        Eigen::Vector2f u1(1.0-Po, Po);
        Eigen::Vector2f u = u1.transpose() *P;
        return u(1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////77
    /////////////////////////////////////////////////////////////////////////////////////////////77
    ///Predict
    /////////////////////////////////////////////////////////////////////////////////////////////77
    /////////////////////////////////////////////////////////////////////////////////////////////77

    float predictOccupancy(int N)
    {
        if( (b_entry_event + b_exit_event)<20.0)
        {
            return 0.5f;
        }
        ///Testing
        float ps = getOccStaticLikelihood();
        if(ps>0.6) return 0.9;
        ///end testing


        float Po = computeShortTermOccupancy();
        float Lex = exitL();
        float Len = entryL();
        Eigen::Matrix2f P;
        P(0,0) = (1.0-Len);
        P(0,1) =  Len;
        P(1,0)	 = Lex;
        P(1,1) = (1-Lex);
        Eigen::Vector2f u1(1.0-Po, Po);
        for(int i=0; i<N; i++) P = P*P;
        Eigen::Vector2f u = u1.transpose() *P;

        return u(1);
    }


};
#endif

