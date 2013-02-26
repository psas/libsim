/**
 * @file
 * @author  John Brewer <http://www.jera.com/techinfo/jtns/jtn002.html>
 *
 * @section LICENSE
 * You may use the code in this tech note for any purpose, with the 
 * understanding that it comes with NO WARRANTY.
 *
 * @brief A (very) minimalist C unit test
 *
 * @section DESCRIPTION
 * A MinUnit test case is just a function that returns 0 (null) if the tests 
 * pass. If the test fails, the function should return a string describing the 
 * failing test. mu_assert is simply a macro that returns a string if the 
 * expression passed to it is false. The mu_runtest macro calls another test 
 * case and returns if that test case fails. That's all there is to it!
 */
#define mu_assert(message, test) do { if (!(test)) return message; } while (0)
#define mu_run_test(test) do { char *message = test(); tests_run++; if (message) return message; } while (0)
extern int tests_run;
