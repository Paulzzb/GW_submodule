% test_log_demo - Demonstrate hierarchical GWlog output with 2-space tags
test()

function test()
QPlog(2, [], 'verbose');       % Verbose level 2 (show all)
QPlog(true, [], 'showtag');% Enable tags

cleanup = QPlog_push('test_log_demo');
QPlog('Starting test_log_demo', 0);

result = func1(3);

QPlog(sprintf('Final result: %f', result), 0);
QPlog('test_log_demo finished', 0);
end


function y = func1(x)
    cleanup = QPlog_push('func1');
    % QPlog('func1', [], '');
    QPlog('Entered func1', 1);

    QPlog('Calling func2...', 1);
    y = func2(x);

    QPlog(sprintf('func2 returned %f', y), 2);
    QPlog('Leaving func1', 1);
end

function out = func2(x)
    % QPlog('func2', [], 'set');
    cleanup = QPlog_push('func2');
    QPlog('Entered func2', 1);

    out = x^2;

    QPlog(sprintf('Result = %f', out), 2);
    QPlog('Leaving func2', 1);
end
