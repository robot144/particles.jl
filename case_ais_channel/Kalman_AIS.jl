#Kalman filter function
function Kalman_filter(i,t,u1,v1)
   R=6371000.0;#Earth radius in meters

   nt=length(t[i].Times)
   
   #Create the time series of the model
   t_model=t[i].Times[1]:Second(1):t[i].Times[nt]
   mt=length(t_model)
   
   dt=1.0;
   count1=1;
   
   Kalman_sec=zeros(length(t[i].Times)-1)
   
   x_k=zeros(6,mt)
   
   C=zeros(2,mt)
   
   #Initial conditions from AIS data
   x_k[1,1]=t[i].Lon[1]
   x_k[2,1]=t[i].Lat[1]
   t_current=(t[i].Times[1]-t0).value/1000
   
   current_x=u1(x_k[1,1],x_k[2,1],0.0,t_current)
   current_y=v1(x_k[1,1],x_k[2,1],0.0,t_current)
   C[1,1]=current_x
   C[2,1]=current_y
   
   #current in direction of heading
   current_along=sind(t[i].HDG[1])*(current_x)+cosd(t[i].HDG[1])*(current_y)                #Calculate the the current along the heading of the ship
   u_ship=(0.514444*(t[i].SOG[1])*(cosd(t[i].COG[1]-t[i].HDG[1])))-current_along
   x_k[3,1]=(u_ship*sind(t[i].HDG[1]))   # u ship clws
   x_k[4,1]=(u_ship*cosd(t[i].HDG[1]))   # v ship clws
   x_k[5,1]=0.0;
   x_k[6,1]=0.0;
   
   
   #Total velocity vectors - saving for output
   stat=zeros(5,mt)
   fill!(stat,NaN)
   utot=zeros(mt)
   vtot=zeros(mt)
   inn_squared=zeros(5,mt)
   fill!(inn_squared,NaN)
   model_stat=zeros(5,mt)
   fill!(model_stat,NaN)
   
   ttrack1=zeros(mt)
   ttrack1[1]=t_current
   utot[1]=current_x+x_k[3,1]+x_k[5,1]
   vtot[1]=current_y+x_k[4,1]+x_k[6,1]
   coslat=zeros(mt)
   fill!(coslat,NaN)
   #Identity matrix
   Id =1.0* Matrix(I,6, 6)
   a=1-(dt/(3600*3));
   sigma_x=sqrt(((1-a^2)*0.3^2)/dt)
   sigma_y=sqrt(((1-a^2)*0.3^2)/dt)
   beta=1;
   #Initial Error Covariance matrix
   P=Diagonal([0.00021^2,0.00013^2,0.0001^2,0.0001^2,0^2,0^2])
   
   #Measurement uncertainty matrix R
   # Lon Lat SOG COG HDG
   Rs= Diagonal([0.00021^2,0.00013^2,(1.7e-2*(180/pi))^2,0.051^2,(1.7e-2*(180/pi))^2])
   
   
   
   
   #Calm water speed of the ship
   for ti=1:mt-1
      #Linearized model x(k+1)=A*x(k)+B*c(k)
      A=[[1 0 (180/(pi*R*cosd(x_k[2,ti])))*dt 0 (180/(pi*R*cosd(x_k[2,ti])))*dt 0]; [0 1 0 (180/(pi*R))*dt 0 (180/(pi*R))*dt];
       [ 0 0 1 0 0 0]; [0 0 0 1 0 0];[0 0 0 0 a 0 ];[0 0 0 0 0 a] ]
      B=[[(180/(pi*R*cosd(x_k[2,ti])))*dt 0];[0  (180/(pi*R))*dt];[0 0]; [0 0];[0 0];[0 0]]
   
       Z=Diagonal([0,0,sqrt(dt),sqrt(dt),sqrt(dt),sqrt(dt)])

      #Forecast
      x_k[:,ti+1]=A*x_k[:,ti]+B*C[:,ti]

       #Noise covariance matrix
       beta1=beta^2*(x_k[3,ti+1]-x_k[3,ti])^2
       beta2=beta^2*(x_k[4,ti+1]-x_k[4,ti])^2
       Qs=Diagonal([0,0,(0.0001^2),(0.0001^2),(sigma_x^2+beta1),(sigma_y^2+beta2)])
   
       P=(A*P*A')+(Z*Qs*Z')


      # Measurement update at ti+1
      time=(t_model[ti+1]-t0).value/1000 #convert to seconds relative to t0
      ttrack1[ti+1]=time
      C[1,ti+1]=u1(x_k[1,ti+1],x_k[2,ti+1],0.0,ttrack1[ti+1])
      C[2,ti+1]=v1(x_k[1,ti+1],x_k[2,ti+1],0.0,ttrack1[ti+1])
      utot[ti+1]=C[1,ti+1]+x_k[3,ti+1]+x_k[5,ti+1]                 # u total velocity (calm water+current+ducurrent)
      vtot[ti+1]=C[2,ti+1]+x_k[4,ti+1]+x_k[6,ti+1]               # v total velocity (calm water+current+dvcurrent)

      # H Matrix
      H=[[1 0 0 0 0 0 ];[0 1 0 0 0 0];
         [0 0 x_k[4,ti+1]/(x_k[3,ti+1]^2+x_k[4,ti+1]^2)  -x_k[3,ti+1]/(x_k[3,ti+1]^2+x_k[4,ti+1]^2) 0 0];
       [0 0 utot[ti+1]/(sqrt(utot[ti+1]^2+vtot[ti+1]^2)) vtot[ti+1]/(sqrt(utot[ti+1]^2+vtot[ti+1]^2)) utot[ti+1]/(sqrt(utot[ti+1]^2+vtot[ti+1]^2)) vtot[ti+1]/(sqrt(utot[ti+1]^2+vtot[ti+1]^2))];
         [0 0 vtot[ti+1]/(utot[ti+1]^2+vtot[ti+1]^2) -utot[ti+1]/(utot[ti+1]^2+vtot[ti+1]^2)  vtot[ti+1]/(utot[ti+1]^2+vtot[ti+1]^2) -utot[ti+1]/(utot[ti+1]^2+vtot[ti+1]^2)]]
   
      #Kalman Gain
      KG=(P*H')*inv((H*P*H'+Rs))
   
      statistics=H*P*H'+Rs
      model_statistics=H*P*H'
       P=0.5*(P+P')

      F=[[x_k[1,ti+1]];
         [x_k[2,ti+1]];
         [atand(x_k[3,ti+1],x_k[4,ti+1])];
         [sqrt(utot[ti+1]^2+vtot[ti+1]^2)];
         [atand(utot[ti+1],vtot[ti+1])]]

      # Analysis-correction
       count2=1
       #Check if the model time matches the observations time
       for count=1:nt
       if (t_model[ti+1]==t[i].Times[count])& (count2==1)
        Y=[[t[i].Lon[count]];
         [t[i].Lat[count]];
         [t[i].HDG[count]];
       [0.51444.*t[i].SOG[count]];
         [t[i].COG[count]]]

          #Create a Date time Series for Kalman (t_model=t_AIS)
         Kalman_sec[count1]=(t_model[ti+1]-t0).value/1000
          #Save the cos of Latitude
         coslat[count1]=cosd(t[i].Lat[count])

         #Save the innovations and the statistics
         innov_squared=(Y-F).^2
               for s=1:5
                   stat[s,ti+1]=statistics[s,s]
                   inn_squared[s,ti+1]=innov_squared[s]
                   model_stat[s,ti+1]=model_statistics[s,s]
               end

        #Convert from degrees to meters
        inn_squared[1,ti+1]=((pi.*R.*coslat[count1])./180).^2*inn_squared[1,ti+1]
        inn_squared[2,ti+1]=((pi.*R)./180).^2*inn_squared[2,ti+1]
        stat[1,ti+1]=((pi.*R.*coslat[count1])./180)^2*stat[1,ti+1]
        stat[2,ti+1]=((pi.*R)./180)^2*stat[2,ti+1]
        model_stat[1,ti+1]=((pi.*R.*coslat[count1])./180)^2*model_stat[1,ti+1]
        model_stat[2,ti+1]=((pi.*R)./180)^2*model_stat[2,ti+1]

       #Update the measurement
      x_k[:,ti+1]=x_k[:,ti+1]+KG*(Y-F) #overwrite x_k

      P=(Id-KG*H)*P

        count1=count1+1
        count2=0
       end
       end
   end #ti
   
   return (x_k,model_stat,stat,inn_squared,utot,vtot,t_model,Kalman_sec,C)
end #Kalman_filter
