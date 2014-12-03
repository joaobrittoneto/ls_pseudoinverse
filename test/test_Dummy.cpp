#include <boost/test/unit_test.hpp>
#include <ls_pseudoinverse/Dummy.hpp>

using namespace ls_pseudoinverse;

BOOST_AUTO_TEST_CASE(it_should_not_crash_when_welcome_is_called)
{
    ls_pseudoinverse::DummyClass dummy;
    dummy.welcome();
}
