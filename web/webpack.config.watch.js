const pkg = require('./package.json');

var path = require('path'),
    webpack = require('webpack'),
    loaders = require('./node_modules/vtk.js/Utilities/config/webpack.loaders.js'),
    plugins = [];
if(process.env.NODE_ENV === 'production') {
    console.log('==> Production build');
    plugins.push(new webpack.DefinePlugin({
        "process.env": {
            NODE_ENV: JSON.stringify("production"),
        },
    }));
}

module.exports = {
  watch: true,
  plugins: plugins,
  entry: pkg.paths.src.js + '/app.js',
  output: {
    // path: pkg.paths.dist.js, // does not work with webpack-stream
    filename: 'bundle.js',
  },
  module: {
      loaders: [
          { test: require.resolve("./src/js/app.js"), loader: "expose-loader?MyWebApp" },
      ].concat(loaders),
    },
    postcss: [
      require('autoprefixer')({ browsers: ['last 2 versions'] }),
    ],
};
