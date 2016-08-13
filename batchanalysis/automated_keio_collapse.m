parent='/home/ashiver/Documents/EMAP/chemical_screen_data/data_camera_iris/';
home={'KEIO1/','KEIO2/','KEIO3/','KEIO4/','KEIO5/','KEIO6/'};
data='0.9.4_iris/';
condition1={'k1b1_new.csv','k2b1_new.csv','k3b1_new.csv','k4b1_new.csv','k5b1_new.csv','k6b1_new.csv'};
condition2={'k1b4_new.csv','k2b4_new.csv','k3b4_new.csv','k4b4_new.csv','k5b4_new.csv','k6b4_new.csv'};
array={'KEIO1_KEY.csv','KEIO2_KEY.csv','KEIO3_KEY.csv','KEIO4_KEY.csv','KEIO5_KEY.csv','KEIO6_KEY.csv'};
structurekey={'a1','b1','a2','b2','a3','b3','a4','b4','a5','b5','a6','b6'};
    
for keio = 1 : 5
    key=structurekey{2*(keio-1)+1};
    [score.(key),datametaC.(key),datamatC.(key),datameta.(key),datamatV.(key)...
        datamatN.(key),datamatS.(key),datamatF.(key),datamat.(key),...
        model.(key),settings.(key)]=scorePHEN([parent,home{keio},condition1{keio}],...
                                [parent,home{keio},array{keio}],...
                                [parent,home{keio},data],...
                                'iris_v0_ecoopacity',{},[],'keio',[]);
    key=structurekey{2*(keio-1)+2};
    [score.(key),datametaC.(key),datamatC.(key),datameta.(key),datamatV.(key)...
        datamatN.(key),datamatS.(key),datamatF.(key),datamat.(key),...
        model.(key),settings.(key)]=scorePHEN([parent,home{keio},condition2{keio}],...
                                [parent,home{keio},array{keio}],...
                                [parent,home{keio},data],...
                                'iris_v0_ecoopacity',{},[],'keio',[]);
end

keio=6;
key=structurekey{11};
[score.(key),datametaC.(key),datamatC.(key),datameta.(key),datamatV.(key)...
        datamatN.(key),datamatS.(key),datamatF.(key),datamat.(key),...
        model.(key),settings.(key)]=scorePHEN([parent,home{keio},condition1{keio}],...
                                [parent,home{keio},array{keio}],...
                                [parent,home{keio},data],...
                                'iris_v0_ecoopacity',{},[],'keio6',keio6ind);
key=structurekey{12};
[score.(key),datametaC.(key),datamatC.(key),datameta.(key),datamatV.(key)...
        datamatN.(key),datamatS.(key),datamatF.(key),datamat.(key),...
        model.(key),settings.(key)]=scorePHEN([parent,home{keio},condition2{keio}],...
                                [parent,home{keio},array{keio}],...
                                [parent,home{keio},data],...
                                'iris_v0_ecoopacity',{},[],'keio6',keio6ind);