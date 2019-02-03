

imageFile = 'NmNf_80';

images_object = matfile(imageFile);%order is: all faces for first frame... all faces for last frame 
mask = images_object.mask;

sdeve = 5.521;
c = 2;
% cs = .0015;
js = 0;

trials = 30000;

cs_values = .0011:.0001:.002;

a = images_object.a;

pc_IO_gender.cs_values = cs_values;

for cs_num = 1:length(cs_values)
    
    cs = cs_values(cs_num);
    pc_IO_gender.(['cs_', num2str(cs_num)]) = IO_gender(a,mask,c,cs,sdeve,js,trials, ...
        'saveDir', 'SavedTemps_hn_Mparams', 'printText', ['cs_num: ', num2str(cs_num)]);
    save('pc_IO_gender', 'pc_IO_gender')
end


    

