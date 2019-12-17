%-------------------------------------------------------
function H = JCBB_test (prediction, observations, compatibility)
% 
%-------------------------------------------------------

Best = zeros(1, observations.m);

Best = JCBB_R (prediction, observations, compatibility, [], 1, Best);

H = Best;
end

%-------------------------------------------------------
function Best = JCBB_R (prediction, observations, compatibility, H, i, Best)
%
% disp('entering function');
if pairings(Best)<100
%-------------------------------------------------------
if i > observations.m % leaf node?
    if pairings(H) > pairings(Best) % did better?
        Best = H;
    end
else
    individually_compatible = find(compatibility.ic(i,:));
    individually_compatible(ismember(individually_compatible,H))=[];
    for j = individually_compatible
        if jointly_compatible(prediction, observations, [H j])
            Best = JCBB_R(prediction, observations, compatibility, [H j], i + 1, Best); %pairing (Ei, Fj) accepted 
        end
    end
%    if pairings(H) + length(compatibility.candidates.observations) - i >= pairings(Best) % can do better?
%    if pairings(H) + observations.m - i >= pairings(Best) % can do better?
    if pairings(H) + pairings(compatibility.AL(i+1:end)) >= pairings(Best) % can do better?
        Best = JCBB_R(prediction, observations, compatibility, [H 0], i + 1, Best); % star node: Ei not paired
    end
end

end

end

%-------------------------------------------------------
% 
%-------------------------------------------------------
function p = pairings(H)

p = length(find(H));
end