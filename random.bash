#random bash script, growing repository of commands 
#for loop for counting. 
file=test.txt
for var in {1..1000000}
do  
echo "$var"
echo "$var" >> $file
done
