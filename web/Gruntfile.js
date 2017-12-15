module.exports = function(grunt) {
  grunt.initConfig({
    css_selectors: {
      options: {
        mutations: [
          {prefix: '#doxygen-content'}
        ]
      },
      your_target: {
        files: {
          'public/docs/doxygen/doxygen-prefixed.css': ['public/docs/doxygen/doxygen.css'],
        },
      },
    },
  });

  grunt.loadNpmTasks('grunt-css-selectors');

  grunt.registerTask('default', ['css_selectors']);
};
