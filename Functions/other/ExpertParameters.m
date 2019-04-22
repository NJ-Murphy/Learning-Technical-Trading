function expertparams = ExpertParameters(C,ell,k,strats)
    W = size(strats,1);
    L = length(ell);
    K = length(k);
    
    %initialise indices for experts
    hlktind1 = 0;  
    hltind = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ---> Going to need to edit this section to learn
    %%% strategies. Instead of using the matching algorithm to
    %%% create the agents

    % strategy loop
    for c = 1:C % clusters (sectors) of stocks 
        for w0 = 1:W % strategies               
            for ell0 = 1:L  %parameter 1 -> n/n1
                ell1 = ell(ell0);   % ----> p.ell will be a vector of chosen lookback periods to try for n/n1 for each strategy

                for k0 = 1:K   %parameter 2 -> n2/lambda  
                    k1 = k(k0); % ----> p.k will be a vector of chosen lookback periods to try for n2 for each strategy

                    %----> Condition to break out of loop if only one param
                    % is required for current strategy
                    if cell2mat(strats(w0,2)) == 0
                        break
                    end                       

                    %----> Continue loop if the k parameter is < or =
                    %to ell since ell==n/n1 cannot be > or = to k
                    if k1 > ell1
                       %Expert index KLrow(w0,k0,ell0)
                       hlktind1 = hlktind1+1;
                       Expind = hlktind1;
                       %Expind = hlktind(hlktind1); 
                       p.indstorekl(Expind) = w0;
                       expertparam2(Expind,:) = [c,w0,ell1,k1];
                    else 
                        continue
                    end 
                end % k

                %---> Compute the agents performance and find the optimal
                % parameters if the strategy only requires one
                % parameter
                if cell2mat(strats(w0,2)) == 0
                   %this will contain the same information as the if
                   %statement in the k loop above
                   hltind = hltind + 1;
                   Expind = hltind;

                   p.indstorel(Expind) = w0; %store the strategy indices
                   % store the expert parameters for the current
                   % expert
                   expertparam1(Expind,:) = [c,w0,ell1,nan];
                end
            end % ell                
        end % w          
    end % c
    
    expertparams = [expertparam2;expertparam1]; 
    
    % EOF
    %----------------------------------------------------------------------