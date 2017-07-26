const path = require('path');
const webpack = require('webpack');
const merge = require('webpack-merge');
const parts = require('./webpack.parts');
const pkg = require('./package.json');

const commonConfig = merge([
  {
    entry: {
      app: pkg.paths.src.js + 'app.js', // PATHS.app,
    },
    output: {
      path: __dirname + pkg.paths.dist.js,
      filename: 'bundle.js'
    },
    module: {
        rules: [
            {
              test: require.resolve("./src/js/app.js"),
              loader: "expose-loader?MyWebApp"
            }
        ].concat(module.exports = parts.loaders()),
      },
    // postcss: [
      // require('autoprefixer')({ browsers: ['last 2 versions'] }),
    // ],
  },
  // parts.lintJavaScript({ include: pkg.paths.src.js + '/app.js' }),
]);

const productionConfig = merge([
  parts.productionEnv(),
  parts.minifyJavaScript({ useSourceMap: true })
]);

const developmentConfig = merge([
  {
    // plugins: [
      // new webpack.NamedModulesPlugin(),
    // ],
  }
]);

module.exports = function(env) {
  if (env === 'production') {
    return merge(commonConfig, productionConfig);
  }

  return merge(commonConfig, developmentConfig);
};
