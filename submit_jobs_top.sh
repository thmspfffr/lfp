#!/bin/sh

# DEFINE THE NUMBER OF JOBS THAT SHOULD BE SUBMITTED
# -------------------
JOBS_IN_PARALLEL=50
# SLEEP_TIME
# TIME_BETWEEN_JOBS
# -------------------

# Submits n jobs to the torque queing system
OUT=$(qstat | grep -E "tpfeffer.* R" -wc)
QUEUED=$(qstat | grep -E "tpfeffer.* Q" -wc)
OUT=$((OUT + QUEUED))
echo -n "How many jobs do you want to submit? "
read NJOBS
echo "You currently have $OUT job(s) running and/or queued..."
echo -n "(s)ubmit all jobs anyway or (q)ueue? "
read choice

case $choice in
  s)
  echo "Submitting jobs..."
  for i in $( seq 1 $NJOBS); do
    let var1=10*$i;
    echo 'Start Job ' $i 'wait for: ' $var1 's'
    qsub -v var="$var1" submit_jobs.sh
  done
  ;;
  q)  
  AVAIL_SLOTS=$(($JOBS_IN_PARALLEL - $OUT))

  # IF THE NUMBER OF AVAILABLE SLOTS IS 0
  # --------------------------------------
  if (($AVAIL_SLOTS < 1)); then
    JOBS_LEFT=$NJOBS

    echo "Submitting 0 jobs immediately..."

    while true; do
      if (($JOBS_LEFT < 1)); then
        break;
      fi
      echo "Waiting for available slots..."   
      sleep 20s
      OUT=$(qstat | grep -E "tpfeffer.* R" -wc)
      QUEUED=$(qstat | grep -E "tpfeffer.* Q" -wc)
      OUT=$((OUT+QUEUED))

      AVAIL_SLOTS=$(($JOBS_IN_PARALLEL - $OUT))
      
      if ((AVAIL_SLOTS > 0)); then
        echo "Starting $AVAIL_SLOTS jobs..."
        for i in $( seq 1 $AVAIL_SLOTS); do
          let var1=10*$i;  
         	qsub -v var="$var1" submit_jobs.sh
        done
        JOBS_LEFT=$(($JOBS_LEFT-$AVAIL_SLOTS))
       fi
    done
    
  # IF THE NUMBER OF AVAILABLE SLOTS IS > 0
  # --------------------------------------
  else
    echo 'Starting $AVAIL_SLOTS jobs...'
    for i in $( seq 1 $AVAIL_SLOTS); do
      let var1=10*$i;
      #echo 'Start Job ' $i 'wait for: ' $var1 's'
      qsub -v var="$var1" submit_jobs.sh
    done
    JOBS_LEFT=$(($NJOBS - $AVAIL_SLOTS))
    if (($JOBS_LEFT > 0)); then
      while true; do
        if (($JOBS_LEFT < 1)); then
          break;
        fi
        echo "Waiting for more available slots..."
        sleep 20s
        OUT=$(qstat | grep -E "tpfeffer.* R" -wc)
        QUEUED=$(qstat | grep -E "tpfeffer.* Q" -wc)
        OUT=$((OUT+QUEUED))

        AVAIL_SLOTS=$(($JOBS_IN_PARALLEL - $OUT))

        if ((AVAIL_SLOTS > 0)); then
          echo 'Starting $AVAIL_SLOTS jobs...'
          for i in $( seq 1 $AVAIL_SLOTS); do
            let var1=10*$i;  
            qsub -v var="$var1" submit_jobs.sh
          done
          JOBS_LEFT=$(($JOBS_LEFT-$AVAIL_SLOTS))
        fi
      done
   	else
      echo "All jobs successfully submitted!"   
    fi
  fi    
    echo "All jobs successfully submitted!"   
  ;;

esac
exit;

  