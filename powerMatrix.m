function P = powerMatrix(PhiQp,PhiQn,PhiFp,PhiFn,w)

    P = (1i*w/2)*([PhiQp'*PhiFp PhiQp'*PhiFn;
        PhiQn'*PhiFp PhiQn'*PhiFn]-...
        [PhiFp'*PhiQp PhiFp'*PhiQn;
        PhiFn'*PhiQp PhiFn'*PhiQn]);

%     P = (1i*w/2)*[PsiQ'*PhiQp'*PsiF*PhiFp PsiQ'*PhiQp'*PsiF*PhiFn;
%         PsiQ'*PhiQn'*PsiF*PhiFp PsiQ'*PhiQn'*PsiF*PhiFn]-...
%         [PsiF'*PhiFp'*PsiQ*PhiQp PsiF'*PhiFp'*PsiQ*PhiQn;
%         PsiF'*PhiFn'*PsiQ*PhiQp PsiF'*PhiFn'*PsiQ*PhiQn];