function test_hybrid_auto
C = hybrid_fixture('HMXL',1,1,3,3,2,1,[2 2 4 7]);
for change = 0:2
    a = C.Args;
    if change == 1, a{10}(1) = a{10}(1)+.2; end
    if change == 2, a{5}(:,2) = a{5}(:,2)+[.1;-.2;.3;.1]; end
    a{9}.GPU = 'cpu'; [fc,gc] = hybrid_likelihood('HMXL',a,true);
    for mode = {'auto','gpu'}
        a{9}.GPU = mode{1};
        for repeat = 1:2
            [f,g,backend] = hybrid_likelihood('HMXL',a,true);
            assert(max(abs(f-fc)) < 1e-8 && max(abs(g-gc),[],'all') < 1e-6);
            fprintf('Hybrid routing %s -> %s, input change %d, call %d\n',mode{1},backend,change,repeat);
            a{10}(1) = a{10}(1)+.03;
            cpuArgs = a; cpuArgs{9}.GPU = 'cpu';
            d = hybrid_prepare('HMXL',cpuArgs);
            [fc,gc] = hybrid_evaluate(d,a{10},2,true,0,[]);
        end
        a = C.Args;
        if change == 1, a{10}(1) = a{10}(1)+.2; end
        if change == 2, a{5}(:,2) = a{5}(:,2)+[.1;-.2;.3;.1]; end
        cpuArgs = a; cpuArgs{9}.GPU = 'cpu';
        d = hybrid_prepare('HMXL',cpuArgs); [fc,gc] = hybrid_evaluate(d,a{10},2,true,0,[]);
    end
end
for bad = {'invalid',struct('unexpected',1)}
    a{9}.GPU = bad{1}; failed = false;
    try, hybrid_likelihood('HMXL',a,true); catch ex, failed = strcmp(ex.identifier,'DCE:Hybrid:GPUOption'); end
    assert(failed,'Invalid GPU options must be rejected.');
end
fprintf('HYBRID_ROUTING_PASSED\n');
end
