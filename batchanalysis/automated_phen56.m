parent='/home/ashiver/Documents/EMAP/chemical_screen_data/data_phenomics_iris/';
home={'phen1/','phen2/','phen3/','phen4/','phen5/','phen6/'};
data='iris/';
condition1={'p1_4120shiv.csv','p2_4120shiv.csv','p3_4120shiv.csv','p4_4120shiv.csv','p5_4120shiv.csv','p6_4120shiv.csv'};
array={'KEIO1_KEY.csv','KEIO2_KEY.csv','KEIO3_KEY.csv','KEIO4_KEY.csv','KEIO5_KEY.csv','KEIO6_KEY.csv'};
structurekey={'p1','p2','p3','p4','p5','p6'};
    
for keio = 5 : 5 
    key=structurekey{keio};
    [scoremat.(key),datametaC.(key),datamatC.(key),datameta.(key),datamatV.(key),...
        datamatN.(key),datamatS.(key),datamatF.(key),datamat.(key),...
        model.(key),settings.(key)]=scorePHEN([parent,home{keio},condition1{keio}],...
                                [parent,home{keio},array{keio}],...
                                [parent,home{keio},data],...
                                'iris_v0_ecoopacity',{},[],'keio',[]);
end

keio=6;
key=structurekey{6};
[scoremat.(key),datametaC.(key),datamatC.(key),datameta.(key),datamatV.(key),...
        datamatN.(key),datamatS.(key),datamatF.(key),datamat.(key),...
        model.(key),settings.(key)]=scorePHEN([parent,home{keio},condition1{keio}],...
                                [parent,home{keio},array{keio}],...
                                [parent,home{keio},data],...
                               'iris_v0_ecoopacity',{},[],'keio6',keio6ind);
save('/home/ashiver/Desktop/phenU.mat');
exit;