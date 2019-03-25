function genCorr = CalcGeneralCorr(jointAB,a)
    a = Columnize(a);
    
    pA = sum(jointAB,2);
    pB = sum(jointAB,1)';
    pAGivenB = bsxfun(@rdivide,jointAB,pB');
    
    varA_pAGivenB = nansum(bsxfun(@times,pAGivenB,a.^2),1)'-nansum(bsxfun(@times,pAGivenB,a),1)'.^2;
    varA_pA = nansum(pA.*a.^2)-sum(pA.*a).^2;
    
    genCorr = 1-sum(varA_pAGivenB.*pB)/varA_pA;
end