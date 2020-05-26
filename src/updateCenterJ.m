function updated_Centers=updateCenterJ(~,normalxi1,Centers,Ni,~,~,~)
Tin=1;
updated_Centers=zeros(Ni,1);
gama=1/16;
    for i=1:Ni
        if Centers(i)>normalxi1(i)        
  
             dcx= gama*((Centers(i)-normalxi1(i)));
             if (Centers(i)-normalxi1(i))>0.1
                 dcx= gama*0.1;
             end
             
             updated_Centers(i,1)=Centers(i)-dcx;
            
            if updated_Centers(i,1)<0
                updated_Centers(i,1)=0;
            end    
        else

            dc= gama*((normalxi1(i)-Centers(i)));
            if (normalxi1(i)-Centers(i))>0.1     %%?0.1
                dcx= gama*0.1;
            end


            updated_Centers(i,1)=Centers(i)+dc;
            
            if updated_Centers(i,1)>Tin
                updated_Centers(i,1)=Tin;
            end    

        end
    end    
    

end