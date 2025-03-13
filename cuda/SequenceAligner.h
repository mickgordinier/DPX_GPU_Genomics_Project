#include <string>

// Base class from which aligners are derived

class SequenceAligner {
    private:
    protected: // <-- makes it so that derived classes can access, but main() cannot
    std::string reference_str;
    std::string query_str;

    public:
    SequenceAligner(const std::string input_reference, const std::string input_query) : reference_str(input_reference), query_str(input_query){}

    // Virtual functions that will be implemented by derived classes
    virtual void init_matrix() = 0;

    virtual void print_matrix() = 0;
    
    virtual void score_matrix() = 0;
    
    virtual void backtrack() = 0;

    virtual void align() = 0;
    
    virtual void print_results() = 0;
};