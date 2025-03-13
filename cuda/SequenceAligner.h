#include <string>

class SequenceAligner {
    private:
    protected:
    std::string reference_str;
    std::string query_str;

    public:
    SequenceAligner(const std::string input_reference, const std::string input_query) : reference_str(input_reference), query_str(input_query){}

    virtual void init_matrix() = 0;

    virtual void print_matrix() = 0;
    
    virtual void score_matrix() = 0;
    
    virtual void backtrack() = 0;

    virtual void align() = 0;
    
    virtual void print_results() = 0;
};