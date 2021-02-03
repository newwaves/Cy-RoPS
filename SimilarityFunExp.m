%%
function Prob = SimilarityFunExp(FeatVA,FeatVB)
% FeatVA : input 1 X 3 vector or N X 3 vector equal to FeatVB
Prob = [];
if size(FeatVA,1) == 1
    NormFV = norm(FeatVA);
    % StdFV = std(FeatVA);
    for i = 1 : 1 : size(FeatVB,1)
        Prob(end+1,:) = exp(- norm(FeatVB(i,:) - FeatVA) );
    end
else
    if size(FeatVA,1) ~= size(FeatVB,1)
        return;
    end
	for i = 1 : 1 : size(FeatVB,1)
        Prob(end+1,:) = exp(- norm(FeatVB(i,:) - FeatVA(i,:))/1 );
    end
end
end