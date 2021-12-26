%% Function of the "Non Uniform Shadowing"

% We constructed the "non uniform shadowing", dividing the square area in
% four parts where the shadowing inside each of them is randomically
% calculated in a certain (determined) range

function sh=shadowing_non_uniforme(x_ue_pos,y_ue_pos,std_db,shadow)

    rng('default')
    for shd=1:numel(shadow)
        [in,on] = inpolygon(x_ue_pos,y_ue_pos,shadow(shd).x,shadow(shd).y);
        if in==1
            switch shd
                case 1
                    sh= std_db + rand*(7-std_db);
                    return
                case 2
                    sh= std_db + rand*(8-std_db);
                    return
                case 3
                    sh= std_db + rand*(10-std_db);
                    return
                case 4
                    sh= std_db + rand*(12-std_db);
                    return
            end
        end
        
        if on==1 
            % if the user is located between two shadowing squares,  
            % we assign a random value
            sh= std_db + rand*(10-std_db);
        end       
    end            
end
