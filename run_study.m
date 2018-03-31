okOrder = [1 2];
okuse_eye = [0 1];

subjectID = input('Enter Subject ID: ','s');
while isempty(subjectID)
    disp('ERROR: no value entered. Please try again.');
    subjectID = input('Enter Subject ID:','s');
end

exp_init = input('Enter Experimenter initials: ','s');
while isempty(exp_init)
    disp('ERROR: no value entered. Please try again.');
    exp_init = input('Enter Experimenter initials:','s');
end

test_comp = input('Enter Computer name: ','s');
while isempty(test_comp)
    disp('ERROR: no value entered. Please try again.');
    test_comp = input('Enter Computer name:','s');
end

order = input('Enter order number (1,2 ; this should be counterbalanced across subjects): ');
while isempty(order) || sum(okOrder==order)~=1
    disp('ERROR: input must be 1 or 2 . Please try again.');
    order = input('Enter order number (1,2 ; this should be counterbalanced across subjects): ');
end

use_eye = input('Are you using the eyetracker? (1=YES, 0=NO): ');
while isempty(use_eye) || sum(okuse_eye==use_eye)~=1
    disp('ERROR: input must be 1 or 0 . Please try again.');
    use_eye = input('Are you using the eyetracker? (1=YES, 0=NO): ');
end

%%MODIFY IF NEEDED
%2 runs of food_rating
for run=1:2
     food_rating(subjectID,run,use_eye);
 end
 %do all the sorting and forming of choice pairs
 sort_ratings(subjectID,order);
 cat_form_probe_pairs(subjectID, order, 2); %2 repetitions of each unique choice pair for CAT_probe
 form_food_choice_pairs(subjectID, order,3); %split 210 choice trials into 3 runs
 
 for block=1:6
      cat_training(subjectID,order,use_eye,block);
 end
   
 ColorDots_practice(subjectID,test_comp,exp_init,use_eye,0,order);
 
cat_probe(subjectID, order, 1, use_eye);

 for run=1:3
     food_choice(subjectID, run, use_eye)
 end
 
 for run=1:3
    ColorDots_test(subjectID,test_comp,exp_init,use_eye,0,run,order)
 end

