%%
function Prob = SimilarityFunCos(FeatVA,FeatVB)
% FeatVA : input 1 X 3 vector or N X 3 vector equal to FeatVB
Prob = [];
if size(FeatVA,1) == 1
    NormFV = norm(FeatVA);
    for i = 1 : 1 : size(FeatVB,1)
        Prob(end+1,:) = abs(dot(FeatVB(i,:),FeatVA) / (norm(FeatVB(i,:)) * NormFV));
    end
else
    if size(FeatVA,1) ~= size(FeatVB,1)
        return;
    end
    for i = 1 : 1 : size(FeatVB,1)
        Prob(end+1,:) = abs(dot(FeatVB(i,:),FeatVA(i,:)) / (norm(FeatVB(i,:)) * norm(FeatVA(i,:))));
    end
end
end